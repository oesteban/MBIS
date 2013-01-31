// --------------------------------------------------------------------------------------
// File:    MultispectralEMFilter.txx
// Date:    Sep 28, 2011
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

#ifndef MULTISPECTRALEMFILTER_TXX_
#define MULTISPECTRALEMFILTER_TXX_


#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>


namespace itk {

/**
 *
 */
template <class TInputComponent, class TProbabilityPixelType>
MultispectralEMFilter< TInputComponent, TProbabilityPixelType >
::MultispectralEMFilter():
 m_MaximumIteration(10),
 m_NumberOfPureTissues(0),
 m_NumberOfClasses(0),
 m_NumberOfComponents(0),
 m_UseMaskedMode(true),
 m_UsePriorProbabilityImages(false),
 m_IsInitialized(false),
 m_UseOnlyClassify(false),
 m_UseExplicitPVModel( false )
{
	// At least 1 inputs is necessary for a vector image.
	this->SetNumberOfRequiredInputs(0);
}

/**
 *
 */
template <class TInputComponent, class TProbabilityPixelType>
MultispectralEMFilter< TInputComponent, TProbabilityPixelType >
::~MultispectralEMFilter()
{
}


template <class TInputComponent, class TProbabilityPixelType>
void MultispectralEMFilter< TInputComponent, TProbabilityPixelType >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "ImageVector:" << std::endl;
  m_Input->Print( os, indent );
}

template <class TInputComponent, class TProbabilityPixelType>
void MultispectralEMFilter< TInputComponent, TProbabilityPixelType >
::Initialize() {

	if ( m_IsInitialized ) { return; }

	// Initialize Number of classes and pure tissues
	m_NumberOfPureTissues = m_NumberOfClasses;
	m_Posteriori.resize(m_NumberOfPureTissues);

	// Initialize sample
	m_Input = this->GetInputVectorImage();
    m_NumberOfComponents = m_Input->GetNumberOfComponentsPerPixel();
	m_PriorIndexOffset = (unsigned int) m_MaskImage.IsNotNull();

	m_Sample = InputSampleType::New();
	m_Sample->SetImage( m_Input );

	// Prepare mask for filters
	if( m_MaskImage.IsNotNull() ) {
		m_Sample->SetMask( m_MaskImage );
	}

	m_Sample->Initialize();

    // Class statistics initialization
	this->InitializeParameters();

	if ( m_UseExplicitPVModel ) {      // If model is B-GHMRF+PV (Bach et al. 2005)
		this->ExtendToExplicitPVModel(); // TODO check what happens with weights if extend to PV model!
	}

	std::cout << "\t* Initialization Parameters: " << std::endl;
	for (size_t i = 0; i<m_InitialParameters.size(); i++)
		std::cout << "\t\t" << m_InitialParameters[i] << std::endl;


	// TODO CHECK sum priors probability = 1.0
	// Check number of priors
	//if ( m_UsePriorProbabilityImages) {
	//	for( unsigned int i = 0; i< m_Priors.size(); i++) {
	//		if ( m_Priors[i].IsNull() )
	//			m_Priors[i] = m_MaskImage;
	//	}
	//}

	// Components initialization
	double uniform = 1.0 / m_NumberOfClasses;
    for ( unsigned int i = 0; i < m_NumberOfClasses ; i++ ) {
    	m_CVector.push_back( ComponentType::New() );
    	(m_CVector[i])->SetSample( m_Sample );

    	typename InputSampleType::ConstIterator end = m_Sample->End();


    	size_t pos = 0;

    	for (typename InputSampleType::ConstIterator it = m_Sample->Begin(); it != end; ++it, pos++) {
    		double weight = (!m_UsePriorProbabilityImages || m_UseExplicitPVModel)?
					uniform:*(m_Priors[i]->GetBufferPointer()+ it.GetInstanceIdentifier() );
    		(m_CVector[i])->SetWeight( pos, weight );
    	}
	}

	// TODO Histogram computation. Find code for this on rev. 278

	// Initial Proportions initialization
	m_InitialProportions.SetSize(m_NumberOfClasses);
	double totalProportion = 0.0;
	for (int i = m_NumberOfClasses - 1; i >= 0; --i) {
		m_InitialProportions[i] = (i > 0) ?
						1.0 / (double) m_NumberOfClasses : 1.0 - totalProportion;
		totalProportion += m_InitialProportions[i];
	}

	// Estimator set-up
	m_Estimator = EstimatorType::New();
	m_Estimator->SetInput(m_Input);
	m_Estimator->SetMaskImage( m_MaskImage );
	m_Estimator->SetMaximumIteration(m_MaximumIteration);
	m_Estimator->SetUseBiasCorrection(m_UseBiasCorrection);
	m_Estimator->SetInitialProportions( m_InitialProportions );
    for ( unsigned int i = 0 ; i < m_NumberOfClasses ; i++ ) {
      (m_CVector[i])->SetParameters( m_InitialParameters[i] );
	  m_Estimator->AddComponent( (typename ComponentType::Superclass*) (m_CVector[i]) );
	}

    m_IsInitialized = true;
}


template <class TInputComponent, class TProbabilityPixelType>
void MultispectralEMFilter< TInputComponent, TProbabilityPixelType >
::GenerateData()
{
	Superclass::GenerateData();

	// Initialize components and parameters
	this->Initialize();

	if( !m_UseOnlyClassify ) {
		std::cout << "\t* Starting EM estimator..." << std::endl;
		try {
			m_Estimator->Update();
		} catch (itk::ExceptionObject & ex ) {
			std::cerr << "Exception caught updating EM-GMM Estimator:" << ex << std::endl;
		}

		std::cout << "\t* EM estimator finished." << std::endl;

		for (unsigned int i=0; i<m_NumberOfClasses; i++) {
			m_InitialParameters[i] = (m_CVector[i])->GetFullParameters();
		}

		//typename InputVectorImageType::ConstPointer bias = m_Estimator->GetCurrentBias();



	}

    // Classification and output generation
    std::cout << "\t* Classification starting..." << std::endl;
    this->Classify();
}

template <class TInputComponent, class TProbabilityPixelType>
void MultispectralEMFilter< TInputComponent, TProbabilityPixelType >
::Classify() {
    typename DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();
    typename ImageClassifierType::Pointer classifier = ImageClassifierType::New();
    classifier->SetDecisionRule( (typename ImageClassifierType::DecisionRuleType*) decisionRule.GetPointer() );
    classifier->SetImage( m_Input );
    classifier->SetNumberOfClasses( m_NumberOfClasses );

    typename ImageClassLabelVectorObjectType::Pointer  classLabelsObject = ImageClassLabelVectorObjectType::New();
    classifier->SetClassLabels( classLabelsObject );

    ImageClassLabelVectorType &  classLabelsVector = classLabelsObject->Get();

    for ( unsigned int i=0; i<m_NumberOfClasses; i++) {
    	unsigned char label = m_PriorIndexOffset + i;
    	classLabelsVector.push_back( label );
    }

    std::cout << "\t* Setting membership functions" << std::endl;
    typename ImageMembershipFunctionVectorObjectType::Pointer membershipFunctionsObject = ImageMembershipFunctionVectorObjectType::New();
        classifier->SetMembershipFunctions( membershipFunctionsObject );

    ImageMembershipFunctionVectorType &  membershipFunctionsVector = membershipFunctionsObject->Get();

    for (unsigned int i = 0; i < m_NumberOfClasses; i++ ){
    	membershipFunctionsVector.push_back( (m_CVector[i])->GetMembershipFunction() );
    }

    classifier->SetMembershipFunctions( membershipFunctionsObject );
    classifier->Update();

    std::cout << "\t* Classification finished, generating output..." << std::endl;

    OutputImagePointer out = classifier->GetOutput();

    // Generate empty data structure for output image
    OutputImagePointer outputPtr = this->GetOutput(0);
    outputPtr->SetRegions( m_Input->GetLargestPossibleRegion() );
    outputPtr->SetSpacing( m_Input->GetSpacing() );
    outputPtr->SetDirection( m_Input->GetDirection() );
    outputPtr->SetOrigin( m_Input->GetOrigin() );
    outputPtr->Allocate();
    outputPtr->FillBuffer( 0.0 );

    // Iterator to copy to output
    OutputImageIterator im_it( outputPtr, outputPtr->GetLargestPossibleRegion());
    im_it.GoToBegin();

    // Generate (copy data to output)
    while ( !im_it.IsAtEnd() ) {
    	if (  m_MaskImage.IsNull() ||
    	     (m_MaskImage.IsNotNull() && m_MaskImage->GetPixel(im_it.GetIndex() )) ) {
    		im_it.Set( out->GetPixel( im_it.GetIndex() ) );
    	}
    	++im_it;
    }

}

template <class TInputComponent, class TProbabilityPixelType>
void MultispectralEMFilter< TInputComponent, TProbabilityPixelType >
::SetPriorProbabilityImage( unsigned int classId,
		    const MultispectralEMFilter< TInputComponent, TProbabilityPixelType >::ProbabilityImageType * image )
{
	if ( classId >= m_NumberOfClasses )
		itkExceptionMacro("Trying to set a prior probability image for a non-existant class (id " << classId << "> lastClassId " << (m_NumberOfClasses-1) << ").");
	m_Priors[classId]= image;
}

template <class TInputComponent, class TProbabilityPixelType>
const typename MultispectralEMFilter< TInputComponent, TProbabilityPixelType >::ProbabilityImageType * MultispectralEMFilter< TInputComponent, TProbabilityPixelType >
::GetPosteriorProbabilityImage( unsigned int label ) {
	if ( label <= 0 || label > m_Posteriori.size() ) {
		itkExceptionMacro("Trying to get a posteriori probability image for a non-existant class (id " << label << "> lastClassId " << (m_NumberOfClasses) << ").");
	}

	this->Initialize();

	unsigned int classId = label -1;

	// Return cached posteriori image
	//if (m_Posteriori[classId].IsNotNull() ) return m_Posteriori[classId];

	m_Estimator->CalculateDensities();
	ProbabilityImagePointer p = ProbabilityImageType::New();

    // Generate empty data structure for output image
    p->SetRegions( this->GetInputVectorImage()->GetLargestPossibleRegion() );
    p->SetSpacing( this->GetInputVectorImage()->GetSpacing() );
    p->SetDirection( this->GetInputVectorImage()->GetDirection() );
    p->SetOrigin( this->GetInputVectorImage()->GetOrigin() );
    p->Allocate();
    p->FillBuffer( 0.0 );

    // Iterator to copy to output
    ProbabilityImageIterator im_it( p, p->GetLargestPossibleRegion());
    im_it.Begin();

    // Generate (copy data to output)
    while ( !im_it.IsAtEnd() ) {
    	ProbabilityPixelType val = m_CVector[classId]->GetWeight( this->GetInputVectorImage()->ComputeOffset(im_it.GetIndex()) );
		if (val > 0.0) {
			im_it.Set(val);
		}
    	++im_it;
    }
    m_Posteriori[classId] = p;

    return p;
}


template <class TInputComponent, class TProbabilityPixelType>
void MultispectralEMFilter< TInputComponent, TProbabilityPixelType >
::SetPriors( const MultispectralEMFilter< TInputComponent, TProbabilityPixelType >::ProbabilityImagesVector& priors ) {
	//if( priors.size() != m_NumberOfPureTissues )
	//	itkExceptionMacro("Number of priors and classes do not match.");

	typename ProbabilityImagesVector::const_iterator it = priors.begin();
	typename ProbabilityImagesVector::const_iterator end = priors.end();

	while( it!= end ) {
		m_Priors.push_back(*it);
		it++;
	}

	m_UsePriorProbabilityImages = true;
}

template <class TInputComponent, class TProbabilityPixelType>
void MultispectralEMFilter< TInputComponent, TProbabilityPixelType >
::InitializeParameters() {
	size_t nElements = (m_NumberOfComponents + 1)* m_NumberOfComponents;
	bool computeCovariance = false;

	// Check number of parameters
	if (m_InitialParameters.size() != m_NumberOfPureTissues ) {
		itkExceptionMacro( "Initial Parameters don't match number of tissue classes");
	}

	// Check number of elements in parameters. Detect missing covariance matrices
	// Exception when array is not consistent
	for ( size_t c = 0; c < m_NumberOfPureTissues; c++ ) {
		unsigned int nComps = m_InitialParameters[c].Size();
		if ( nComps == m_NumberOfComponents ) {
			computeCovariance = true;
			break;
		 } else {
			if ( nComps != nElements )
				itkExceptionMacro( "Initial Parameters (" << nComps << ") don't match number of channels (" << nElements << ")" );
		}
	}

	// Compute spherical approximated covariance (diagonal matrix) if it's missing
	if ( computeCovariance ) {
		ParametersVectorType tmpParam = m_InitialParameters;

		InputComponentPixelType minPixel = NumericTraits< InputComponentPixelType >::min();
		InputComponentPixelType maxPixel = NumericTraits< InputComponentPixelType >::max();

		MeasurementVectorType max( m_NumberOfComponents );
		max.Fill( minPixel );
		MeasurementVectorType min( m_NumberOfComponents );
		min.Fill( maxPixel );

		typename InputSampleType::ConstIterator it = m_Sample->Begin();
		typename InputSampleType::ConstIterator end = m_Sample->End();
		size_t mId = 0;
		while (it!=end) {
			if( it.GetFrequency() > 0.0 ) {
				const MeasurementVectorType val = it.GetMeasurementVector();
				for (size_t c = 0; c < m_NumberOfComponents; c++) {
					if ( val[c] > max[c] ) max[c] = val[c];
					if ( val[c] < min[c] ) min[c] = val[c];
				}
				mId++;
			}
			++it;
		}

		m_InitialParameters.clear();

		for (unsigned int i = 0; i< m_NumberOfPureTissues; i++) {
			ParametersType statsArray( nElements );
			MeasurementVectorType sigma( m_NumberOfComponents );
			itk::VariableSizeMatrix< double > cov( m_NumberOfComponents , m_NumberOfComponents );

			for( unsigned int comp = 0; comp < m_NumberOfComponents; comp++ ) {
				size_t idx = i * (m_NumberOfComponents) + comp ;
				statsArray.SetElement(comp, tmpParam[i][comp] );

				//double max = m_Histogram->Quantile( comp, 0.9 );
				//double min = m_Histogram->Quantile( comp, 0.15 );

				//sigma[comp] = (max[comp] - min[comp]) * 0.015 * m_NumberOfComponents* m_NumberOfComponents;
				sigma[comp] = (max[comp] - min[comp]) * 0.015;
			}

			for (size_t k = 0; k<m_NumberOfComponents; k++) {
				for (size_t v = 0; v<m_NumberOfComponents; v++) {
					cov(k,v) = (k==v)?sigma[k]*sigma[v]:0;
				}
			}

			typedef vnl_diag_matrix<double>::iterator DiagonalIterator;
			typedef vnl_symmetric_eigensystem<double> Eigensystem;
			vnl_matrix< double > S = cov.GetVnlMatrix();
			Eigensystem* e = new Eigensystem( S );

			DiagonalIterator itD = e->D.begin();
			while ( itD!= e->D.end() ) {
				if (*itD < 0) *itD = 0.;
				itD++;
			}

			vnl_matrix<double> newCov = e->recompose();

			//the determinant is then costless this way
			vnl_matrix_inverse< double > inv_cov( newCov );
			double det = inv_cov.determinant_magnitude();

			typedef vnl_matrix<double>::iterator MatrixIterator;
			MatrixIterator itNS = newCov.begin();
			size_t el = m_NumberOfComponents;
			while( itNS!= newCov.end() ) {
				statsArray.SetElement( el, *itNS );
				el++;
				itNS++;
			}
			// Insert pure tissue class initial parameters
			m_InitialParameters.push_back( statsArray );
		}
	}
}

template <class TInputComponent, class TProbabilityPixelType>
void MultispectralEMFilter< TInputComponent, TProbabilityPixelType >
::ExtendToExplicitPVModel() {
	// Update number of classes
	m_NumberOfClasses += vcl_ceil( m_NumberOfClasses/2.0 );
	size_t nElements = (m_NumberOfComponents + 1)* m_NumberOfComponents;

	// Generate and insert PV classes initial parameters
	for ( size_t i = m_NumberOfPureTissues - 1; i > 0; i--) {
		// Create objects to keep the PV initial parameters
		ParametersType PVArray( nElements );

		itk::Array<double> val1 = m_InitialParameters[i-1];
		itk::Array<double> val2 = m_InitialParameters[i];

		for( unsigned int comp = 0; comp < nElements; comp++ ) {
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
#endif /* MULTISPECTRALEMFILTER_TXX_ */
