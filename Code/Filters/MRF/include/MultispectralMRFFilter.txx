// --------------------------------------------------------------------------------------
// File:    MultispectralMRFFilter.txx
// Date:    Nov 9, 2011
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of MBIS
//
// MBIS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MBIS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MBIS.  If not, see <http://www.gnu.org/licenses/>.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef MULTISPECTRALMRFFILTER_TXX_
#define MULTISPECTRALMRFFILTER_TXX_

#include <boost/filesystem.hpp>
#include "itkEventObject.h"
#include "MultispectralMRFFilter.h"

namespace bfs = boost::filesystem;

namespace itk {

/**
 *
 */
template <class TInputComponent, class TProbabilityPixelType, class TEnergyValType>
MultispectralMRFFilter< TInputComponent, TProbabilityPixelType, TEnergyValType >
::MultispectralMRFFilter():
 m_NumberOfClasses(3),
 m_NumberOfPureTissues(3),
 m_Lambda(0.5),
 m_ExternalLambda(1.0),
 m_MRFIterations(1),
 m_CurrentIteration(0),
 m_UseOutputPriors(false),
 m_UseExplicitPVModel(false),
 m_UseCachedEnergyMap(true),
 m_MinimizationMode (GCOptimizerType::EXPANSION),
 m_OutputPrefix("")
{
	this->SetNumberOfRequiredInputs(0);
	m_EnergyFunction = EnergyType::New();
}

/**
 *
 */
template <class TInputComponent, class TProbabilityPixelType, class TEnergyValType>
MultispectralMRFFilter< TInputComponent, TProbabilityPixelType, TEnergyValType >
::~MultispectralMRFFilter()
{
}


template <class TInputComponent, class TProbabilityPixelType, class TEnergyValType>
void MultispectralMRFFilter< TInputComponent, TProbabilityPixelType, TEnergyValType >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "ImageVector:" << std::endl;
  m_Input->Print( os, indent );
}


template <class TInputComponent, class TProbabilityPixelType, class TEnergyValType>
void MultispectralMRFFilter< TInputComponent, TProbabilityPixelType, TEnergyValType >
::GenerateData()
{
	// Initial general settings -----------------------------------------------------------------------------
	m_Input = this->GetInputVectorImage();
    m_NumberOfComponents = m_Input->GetNumberOfComponentsPerPixel();
    m_NumberOfPureTissues = m_NumberOfClasses;
	m_MatrixOffset = (m_MaskImage.IsNotNull())?1:0;
    m_NumberOfClasses = ( !m_UseExplicitPVModel )?m_NumberOfPureTissues:m_NumberOfPureTissues+vcl_ceil( m_NumberOfPureTissues/2.0 );

	// Initialize output ------------------------------------------------------------------------------------
	OutputImagePointer outputPtr = this->GetOutput(0);
    outputPtr->SetRegions( m_Input->GetLargestPossibleRegion() );
    outputPtr->CopyInformation( m_Input );
    outputPtr->Allocate();
    outputPtr->FillBuffer( 0.0 );
    outputPtr->Update();

	// Initialize Sample for non spatially-constrained statistical analysis ---------------------------------
	m_Sample = InputSampleType::New();
	m_Sample->SetImage( m_Input );
	if( m_MaskImage.IsNotNull() ) {
		m_Sample->SetMask( m_MaskImage );
	}
	m_Sample->Initialize();
	size_t sampleSize = m_Sample->Size();

	// Model initialization ---------------------------------------------------------------------------------
	if (m_UseExplicitPVModel) this->ExtendToExplicitPVModel();
	this->InitializeModel();

	// Costs matrix initialization block --------------------------------------------------------------------
	if ( bfs::exists( m_MatrixFile ) ) this->LoadMatrixFile();
	std::cout << "\t* Smoothness energy costs matrix:" << std::endl;
	for ( size_t i = 0; i < m_EnergyMatrix.Rows(); i++  ) {
		std::cout << "\t\t[";
		for ( size_t j = 0; j < m_EnergyMatrix.Cols(); j++  ) {
			if (m_EnergyMatrix(i,j)<1e3)
				std::cout << std::setw(5) << std::setprecision(3) << m_EnergyMatrix(i,j);
			else
				std::cout << std::setw(5) << "inf";
		}
		std::cout << "]" << std::endl;
	}
	// End costs matrix initialization block ----------------------------------------------

	m_EnergyFunction->SetOutput( outputPtr );
	m_EnergyFunction->SetInput( m_Input );
	if (m_MaskImage.IsNotNull() ) m_EnergyFunction->SetMaskImage( m_MaskImage );
	m_EnergyFunction->SetNumberOfClasses( m_NumberOfClasses );
	m_EnergyFunction->SetLambda( m_Lambda );
	m_EnergyFunction->SetExternalLambda( m_ExternalLambda );
	m_EnergyFunction->SetEnergyMatrix( m_EnergyMatrix );

	m_EnergyFunction->Initialize();

	// Set-up external energies field (Vx) --------------------------------------------------
	if ( m_SpatialEnergyMap.size() > 0) {
		this->InitializeExternalFieldMaps();
	}
	// End Set-up external energies field (Vx) ----------------------------------------------

	for(size_t i = 0; i< m_NumberOfClasses; i++) {
		m_EnergyFunction->SetModelParametersForLabel( i+1, m_InitialParameters[i] );
	}

	this->InvokeEvent( StartEvent() );


    while( m_CurrentIteration < m_MRFIterations ) {
    	std::cout << "\t* Iteration " << (int) m_CurrentIteration << "." << std::flush;

    	if( m_CurrentIteration>0 ) { // re-estimate the model;

    		typename ParameterEstimator::Pointer est = ParameterEstimator::New();
    		est->SetInputVectorImage( m_Input );
    		est->SetMaskImage( m_MaskImage );


    		std::cout << "Re-estimated parameters:" << std::endl;
    		for(size_t i = 0; i< m_NumberOfClasses; i++) {
        		ProbabilityImageConstPointer prior = m_EnergyFunction->GetPostProbMapForLabel(i+1);
    			est->SetPrior( prior );
        		est->Update();

        		m_InitialParameters[i] = est->GetOutputParameters();
        		m_EnergyFunction->SetModelParametersForLabel( i+1, m_InitialParameters[i] );
    		}
    	}
    	else {
    		std::cout << "Initial Parameters:" << std::endl;
    	}

    	for( size_t i = 0; i< m_NumberOfClasses; i++) {
    		std::cout << "\t\t" << m_InitialParameters[i] << std::endl;
    	}


    	this->MinimizeFunctional();
    	m_EnergyFunction->ComputePostProbMaps();
    	this->InvokeEvent( IterationEvent() );
    	m_CurrentIteration++;
    }

    this->InvokeEvent( EndEvent() );
}

template <class TInputComponent, class TProbabilityPixelType, class TEnergyValType>
void MultispectralMRFFilter< TInputComponent, TProbabilityPixelType, TEnergyValType >::MinimizeFunctional()
{
	typename OptimizerType::Pointer o;
	switch (m_MinimizationMode ) {
	case GCOptimizerType::SWAP :
		o = SwapOptimizerType::New();
		std::cout << "\t\tMinimization algorithm is a-b-Swap" << std::endl;
		break;
	default:
		o = ExpansionOptimizerType::New();
		std::cout << "\t\tMinimization algorithm is a-Expansion" << std::endl;
		break;
	}
	o->SetMRFEnergy( m_EnergyFunction );
	o->StartOptimization();

}


template <class TInputComponent, class TProbabilityPixelType, class TEnergyValType>
void MultispectralMRFFilter< TInputComponent, TProbabilityPixelType, TEnergyValType >
::InitializeModel() {
	// Components initialization
	double uniform = 1.0 / m_NumberOfClasses;
	size_t sampleSize = m_Sample->Size();

	for ( size_t i = 0; i < m_NumberOfClasses; i++ ) {
		m_CVector.push_back( ComponentType::New() );
		(m_CVector[i])->SetSample( m_Sample );

		for (size_t pos=0; pos < sampleSize; pos++) {
	    	(m_CVector[i])->SetWeight( pos, uniform );
	    }

		(m_CVector[i])->SetParameters( m_InitialParameters[i] );
	}
}

template <class TInputComponent, class TProbabilityPixelType, class TEnergyValType>
void MultispectralMRFFilter< TInputComponent, TProbabilityPixelType, TEnergyValType >
::ExtendToExplicitPVModel() {
	// Update number of classes
	size_t nElements = (m_NumberOfComponents + 1)* m_NumberOfComponents;

	if( m_InitialParameters.size() == m_NumberOfPureTissues ) {

		// Generate and insert PV classes initial parameters
		for ( size_t i = m_NumberOfPureTissues - 1; i > 0; i--) {
			// Create objects to keep the PV initial parameters
			ParametersType PVArray( nElements );

			itk::Array<double> val1 = m_InitialParameters[i-1];
			itk::Array<double> val2 = m_InitialParameters[i];

			for( size_t comp = 0; comp < nElements; comp++ ) {
				// Compute PV means
				double value = ( val1[comp] + val2[comp])*0.5;
				if ( comp >= m_NumberOfComponents )
					value = value * 0.5;
				// Compute PV covariance matrix
				PVArray.SetElement(comp, value);
			}


			typename ParametersVectorType::iterator it = m_InitialParameters.begin();
			std::advance( it, i );
			// Insert PV class initial parameters
			m_InitialParameters.insert( it, PVArray );
		}
	}
}


template <class TInputComponent, class TProbabilityPixelType, class TEnergyValType>
void MultispectralMRFFilter< TInputComponent, TProbabilityPixelType, TEnergyValType >
::LoadMatrixFile( ){
	std::ifstream file(m_MatrixFile.c_str());
	std::string line;
	size_t nClasses = m_NumberOfClasses + m_MatrixOffset;
	m_EnergyMatrix = EnergyMatrixType( nClasses, nClasses );

	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep(",;| []");

	size_t row = 0;
	size_t col = 0;


	while ( getline(file, line)) {
		size_t pos = line.find("#");

		if( pos!=0 ) {
			line = line.substr(0,pos);

			tokenizer toker(line, sep);
			tokenizer::iterator iterator = toker.begin();

			while (iterator != toker.end()) {
				float val = 0.0;
				std::string item = boost::lexical_cast< std::string >( *iterator++ );

				if (!item.compare("inf") ) {
					val = MRF_INF;
				} else {
					val = boost::lexical_cast<float>( item );
				}
				m_EnergyMatrix(row,col) = val;
				col = (col+ 1)%( nClasses );
			}
			row++;
		}
	}

	file.close();
}

template <class TInputComponent, class TProbabilityPixelType, class TEnergyValType>
void MultispectralMRFFilter< TInputComponent, TProbabilityPixelType, TEnergyValType >
::InitializeExternalFieldMaps() {
	ProbabilityImagesVector maps;
	maps.resize( m_NumberOfClasses );

	std::vector< LabelType >::iterator begin = m_SpatialEnergyMapIndex.begin();
	std::vector< LabelType >::iterator end = m_SpatialEnergyMapIndex.end();
	std::vector< LabelType >::iterator el;

	// Initialize non-setted maps
	for( size_t i = 0; i < m_NumberOfClasses; i++ ) {
		el = std::find( begin, end, (LabelType) (i+1) );

		if( el != end ) {
			maps[i] = m_SpatialEnergyMap[el-begin];
		} else {
			maps[i] = ProbabilityImageType::New();
			maps[i]->CopyInformation( m_Input );
			maps[i]->SetRegions( m_Input->GetRequestedRegion() );
			maps[i]->Allocate();
			maps[i]->FillBuffer(1.0);
			maps[i]->Update();
		}
	}

	// Normalize
	size_t last = m_Input->GetLargestPossibleRegion().GetNumberOfPixels();
	ProbabilityImagePointer totalProb = ProbabilityImageType::New();
	totalProb->CopyInformation( m_Input );
	totalProb->SetRegions( m_Input->GetRequestedRegion() );
	totalProb->Allocate();
	totalProb->FillBuffer(0.0);
	totalProb->Update();
	ProbabilityPixelType* totalBuffer = totalProb->GetBufferPointer();

	float unifProb = (1.0/m_NumberOfClasses);

	for( size_t cl = 0; cl < m_NumberOfClasses; cl++ ) {
		size_t pix;
		ProbabilityPixelType* buffer = maps[cl]->GetBufferPointer();

		for( size_t pix = 0; pix < last; pix++) {
			(*(totalBuffer + pix)) += (*(buffer + pix));
		}
	}

	for( size_t cl = 0; cl < m_NumberOfClasses; cl++ ) {
		size_t pix;
		ProbabilityPixelType* buffer = maps[cl]->GetBufferPointer();

		for( size_t pix = 0; pix < last; pix++) {
			ProbabilityPixelType val = *( totalBuffer + pix );
			*(buffer + pix)= ( val > 0 )?
					((*(buffer + pix))/val) : unifProb;
		}
	}

	// Set maps
	for( size_t i = 0; i < m_NumberOfClasses; i++ ) {
		m_EnergyFunction->SetSpatialEnergyMap( i+1 , maps[i] );
	}

}

template <class TInputComponent, class TProbabilityPixelType, class TEnergyValType>
void MultispectralMRFFilter< TInputComponent, TProbabilityPixelType, TEnergyValType >
::WriteIteration( )
 {
	std::stringstream filename;
	filename << m_OutputPrefix <<  "_mrf_it" << m_CurrentIteration << ".nii.gz";
	typename itk::ImageFileWriter<OutputImageType>::Pointer w = itk::ImageFileWriter<OutputImageType>::New();
	w->SetInput( this->GetOutput(0) );
	w->SetFileName( filename.str() );
	w->Update();
 }

template <class TInputComponent, class TProbabilityPixelType, class TEnergyValType>
void MultispectralMRFFilter< TInputComponent, TProbabilityPixelType, TEnergyValType >
::WriteMaps( )
 {
	std::stringstream filename;

	typename itk::ImageFileWriter<ProbabilityImageType>::Pointer w2 = itk::ImageFileWriter<ProbabilityImageType>::New();
	for ( size_t i = 0; i < m_NumberOfClasses; i++) {
		filename.str("");
		filename << m_OutputPrefix << "_it" << m_CurrentIteration <<  "_mrf_pmap" << i << ".nii.gz";

		w2->SetInput( m_EnergyFunction->GetDataProbMapForLabel(i+1) );
		w2->SetFileName( filename.str() );
		w2->Update();

		filename.str("");
		filename << m_OutputPrefix << "_it" << m_CurrentIteration <<  "_mrf_emap" << i << ".nii.gz";

		w2->SetInput( m_EnergyFunction->GetDataEnergyMapForLabel(i+1) );
		w2->SetFileName( filename.str() );
		w2->Update();

		filename.str("");
		filename << m_OutputPrefix << "_it" << m_CurrentIteration <<  "_mrf_ppmap" << i << ".nii.gz";

		w2->SetInput( m_EnergyFunction->GetPostProbMapForLabel(i+1) );
		w2->SetFileName( filename.str() );
		w2->Update();

		filename.str("");
		filename << m_OutputPrefix << "_it" << m_CurrentIteration <<  "_mrf_semap" << i << ".nii.gz";

		w2->SetInput( m_EnergyFunction->GetSpatialEnergyMapForLabel(i+1) );
		w2->SetFileName( filename.str() );
		w2->Update();

		filename.str("");
		filename << m_OutputPrefix << "_it" << m_CurrentIteration <<  "_mrf_MRFpmap" << i << ".nii.gz";

		w2->SetInput( m_EnergyFunction->GetMRFProbMapForLabel(i+1) );
		w2->SetFileName( filename.str() );
		w2->Update();
	}

	filename.str("");
	filename << m_OutputPrefix << "_it" << m_CurrentIteration <<  "_mrf_mrfmap.nii.gz";
	m_EnergyFunction->SaveMRFEnergyMap(filename.str());
 }

}


#endif /* MULTISPECTRALMRFFILTER_TXX_ */
