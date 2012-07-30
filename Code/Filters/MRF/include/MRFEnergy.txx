// --------------------------------------------------------------------------------------
// File:    MRFEnergy.txx
// Date:    Mar 16, 2012
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

#ifndef MRFENERGY_TXX_
#define MRFENERGY_TXX_


#include "MRFEnergy.h"

template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
const float MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
MAX_LOG = - vcl_log( vnl_math::eps );

template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
MRFEnergy() :
   m_NumberOfClasses(0),
   m_Lambda(1.0),
   m_ExternalLambda(1.0),
   m_MatrixOffset(1),
   m_IsInitialized(false),
   m_TotalDataEnergy(0.0),
   m_TotalSmoothnessEnergy(0.0),
   m_OutputBuffer(0),
   m_DataEnergyMapBuffer(0),
   m_ModelModified(true),
   m_LabelingModified(true)
{
	//m_EnergyMatrix = NULL;
};

template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
void MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
Initialize() {
	m_MatrixOffset = (m_MaskImage.IsNotNull())?1:0;

	// Initialize default energy matrix (Pott's Model)
	if( m_EnergyMatrix.Rows() == 0 || m_EnergyMatrix.Cols() == 0 ) {
		m_EnergyMatrix = EnergyMatrixType( m_NumberOfClasses + m_MatrixOffset, m_NumberOfClasses + m_MatrixOffset);
		m_EnergyMatrix.Fill(0.0);

		for (size_t i = m_MatrixOffset; i<m_EnergyMatrix.Cols(); i++) {
			for (size_t j = m_MatrixOffset; j<m_EnergyMatrix.Rows(); j++) {
				m_EnergyMatrix(i,j) = (float) ( i!=j );
			}
		}
	}

	m_EnergyMatrix *= this->m_Lambda;

	m_MembershipVector.resize( m_NumberOfClasses );
	for ( size_t i = 0; i < m_NumberOfClasses; i++)
		m_MembershipVector[i] = DefaultMembershipFunction::New();

	// Initialize Data Energy Maps
	for( size_t i = 0; i<m_NumberOfClasses; i++) {
		m_DataEnergyMap.push_back( EnergyMapType::New() );
		m_DataEnergyMap[i]->SetRegions( m_Output->GetRequestedRegion() );
		m_DataEnergyMap[i]->CopyInformation( m_Output );
		m_DataEnergyMap[i]->Allocate();
		m_DataEnergyMap[i]->FillBuffer(0.0);
		m_DataEnergyMap[i]->Update();
		m_DataEnergyMapBuffer.push_back( m_DataEnergyMap[i]->GetBufferPointer() );

		m_DataProbMap.push_back(  EnergyMapType::New() );
		m_DataProbMap[i]->CopyInformation( m_Output );
		m_DataProbMap[i]->SetRegions( m_Output->GetRequestedRegion() );
		m_DataProbMap[i]->Allocate();
		m_DataProbMap[i]->FillBuffer(0.0);
		m_DataProbMap[i]->Update();
		m_DataProbMapBuffer.push_back( m_DataProbMap[i]->GetBufferPointer() );

		m_PostProbMap.push_back(  EnergyMapType::New() );
		m_PostProbMap[i]->CopyInformation( m_Output );
		m_PostProbMap[i]->SetRegions( m_Output->GetRequestedRegion() );
		m_PostProbMap[i]->Allocate();
		m_PostProbMap[i]->FillBuffer(0.0);
		m_PostProbMap[i]->Update();
		m_PostProbMapBuffer.push_back( m_PostProbMap[i]->GetBufferPointer() );

		m_MRFProbMap.push_back(  EnergyMapType::New() );
		m_MRFProbMap[i]->CopyInformation( m_Output );
		m_MRFProbMap[i]->SetRegions( m_Output->GetRequestedRegion() );
		m_MRFProbMap[i]->Allocate();
		m_MRFProbMap[i]->FillBuffer(0.0);
		m_MRFProbMap[i]->Update();
		m_MRFProbMapBuffer.push_back( m_MRFProbMap[i]->GetBufferPointer() );

		m_SpatialEnergyMap.push_back(  EnergyMapType::New() );
		m_SpatialEnergyMap[i]->CopyInformation( m_Output );
		m_SpatialEnergyMap[i]->SetRegions( m_Output->GetRequestedRegion() );
		m_SpatialEnergyMap[i]->Allocate();
		m_SpatialEnergyMap[i]->FillBuffer(0.0);
		m_SpatialEnergyMap[i]->Update();
		m_SpatialEnergyMapBuffer.push_back( m_SpatialEnergyMap[i]->GetBufferPointer() );
	}
	m_OutputBuffer = m_Output->GetBufferPointer();

	// Initialize Sample for non spatially-constrained statistical analysis ---------------------------------
	m_Sample = InputSampleType::New();
	m_Sample->SetImage( m_Input );
	if( m_MaskImage.IsNotNull() ) {
		m_Sample->SetMask( m_MaskImage );
	}
	m_Sample->Initialize();

	// Initialize m_IdHolder (SitesId holder)
	m_IdHolder = m_Sample->GetIdHolder();


	// Initialize OutputBuffer with label 1
	LabelType* lBuffer = this->m_OutputBuffer;
	typename InputSampleType::ConstIterator end = m_Sample->End();
   	for (typename InputSampleType::ConstIterator it = m_Sample->Begin(); it != end; ++it) {
		*( lBuffer + it.GetInstanceIdentifier() ) = 1.0;
	}


	m_IsInitialized = true;
}

template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
void MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
UpdateTotalSmoothnessEnergy() {
	EnergyValType eng = 0;
	size_t offset;

	typename InputSampleType::ConstIterator end = m_Sample->End();
   	for (typename InputSampleType::ConstIterator it = m_Sample->Begin(); it != end; ++it) {
   		offset = it.GetInstanceIdentifier();
   		eng += this->GetCliqueEnergy( offset );
	}

	this->m_TotalSmoothnessEnergy = eng;
	this->m_LabelingModified = false;
}

template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
void MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
UpdateTotalDataEnergy() {
	EnergyValType eng = 0;
	size_t offset;

	typename InputSampleType::ConstIterator end = m_Sample->End();
   	for (typename InputSampleType::ConstIterator it = m_Sample->Begin(); it != end; ++it) {
   		offset = it.GetInstanceIdentifier();
		eng+= this->GetSampleEnergy( offset );
	}
	this->m_TotalDataEnergy = eng;
}


template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
typename MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::EnergyValType
MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
GetCliqueEnergy(size_t offset) {
	LabelType* lBuffer = this->m_Output->GetBufferPointer();
	LabelType cLabel = *( lBuffer + offset);
	assert( cLabel!= 0 );

	return this->GetCliqueEnergyForLabel(offset, cLabel);
}

template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
typename MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::EnergyValType
MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
GetCliqueEnergyForLabel(size_t offset, LabelType cLabel) {
	assert( cLabel > 0 );
	size_t nClasses = this->m_NumberOfClasses;
	EnergyValType Ux = 0.0;

	typename LabelImageType::IndexType centerIdx = this->m_Output->ComputeIndex( offset );
	typename LabelImageType::SpacingType spacing = this->m_Output->GetSpacing();
	typename LabelImageType::SizeType size = this->m_Output->GetLargestPossibleRegion().GetSize();
	typename LabelImageType::IndexType nIdx;
	LabelType nLabel;
	LabelType* lBuffer = this->m_Output->GetBufferPointer();
	size_t nOffset;

	for (size_t dim = 0; dim < ImageDimension; dim++) {
		for( size_t n = 0; n<2; n++ ) {
			nIdx = centerIdx;
			nIdx[dim] += ( -2*n + 1 );
			if ( nIdx[dim] >= 0 && nIdx[dim] < size[dim] ) {
				nOffset = this->m_Output->ComputeOffset(nIdx);
				nLabel = *( lBuffer + nOffset );
				Ux += ( this->GetNeighborCost(cLabel, nLabel) / spacing[dim] );
			}
		}

	}
	if( Ux > MRF_INF ) Ux = MRF_INF;

	return Ux;
}

template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
inline typename MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::EnergyValType
MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
GetExternalFieldEnergy(size_t offset) {
	LabelType* lBuffer = this->m_Output->GetBufferPointer();
	LabelType cLabel = *( lBuffer + offset);
	assert( cLabel!= 0 );

	return *(m_SpatialEnergyMapBuffer[cLabel-1] + offset);
}

template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
inline typename MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::EnergyValType
MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
GetExternalFieldEnergyForLabel(size_t offset, LabelType cLabel) {
	return *(m_SpatialEnergyMapBuffer[cLabel-1] + offset);
}

template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
void MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
UpdateMaps() {
	if (!m_IsInitialized) Initialize();

	EnergyValType Gt;
	EnergyValType Ed;
	MeasurementVectorType m;
	size_t measurementVectorSize = m_Sample->GetMeasurementVectorSize();
	size_t offset;

	typename InputSampleType::ConstIterator end = m_Sample->End();


	for (size_t labelId = 0; labelId < m_NumberOfClasses; labelId++ ) {
		for (typename InputSampleType::ConstIterator it = m_Sample->Begin(); it != end; ++it) {
			m = it.GetMeasurementVector();
			offset = it.GetInstanceIdentifier();
			Ed = m_MembershipVector[labelId]->EvaluateEnergy(m);
			*(m_DataEnergyMapBuffer[labelId] + offset) = ( Ed < DATA_MAX_ENERGY )?Ed:DATA_MAX_ENERGY;
   		}
	}

   	this->UpdateTotalDataEnergy();

   	m_ModelModified = false;
}

template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
void MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
ComputePostProbMaps() {
	UpdateMaps();

	size_t nClasses = m_NumberOfClasses;
	double uniformP = 1.0 / nClasses;

    // Iterate over every active (inside mask) voxel in the image
	MeasurementVectorType m;
	size_t offset;

	EnergyValType Zx, GtTotal, postPTotal, normP, sumVij, Vx, Ux;
	itk::Array<EnergyValType> Px, Gt, postP;
	Px.SetSize(nClasses);
	Gt.SetSize(nClasses);
	postP.SetSize(nClasses);

	typename InputSampleType::ConstIterator end = m_Sample->End();
   	for (typename InputSampleType::ConstIterator it = m_Sample->Begin(); it != end; ++it) {
   		offset = it.GetInstanceIdentifier();
   		m = it.GetMeasurementVector();
		Zx = 0.0;
		GtTotal = 0.0;
		postPTotal= 0.0;

   		// Perform mrf prior calculation of pixel
   		for ( size_t c = 0; c < nClasses; c++ ){
   			// Gt
   			Gt[c] =  m_MembershipVector[c]->Evaluate(m);
   			GtTotal += Gt[c];

   			// Ptlx (Bach 2005, form. 15 )
   			Vx     = this->GetExternalFieldEnergyForLabel( offset, c+1 );
   			sumVij = GetCliqueEnergyForLabel(offset, c+1);
   			*(m_MRFProbMapBuffer[c] + offset ) = sumVij;
   			Ux = Vx + sumVij;
   			if ( Ux > MRF_INF ) Ux = MRF_INF;
   			Px[c] = vcl_exp( (-Ux) );

   			// Find Z in this pixel for normalization (Bach 2005, form. 19)
   			Zx += Px[c];
   		}

   		for ( size_t c = 0; c < nClasses; c++ ){
			// Ptly * Pty = Gt * Ptlx
   			double gt = uniformP;
   			if( GtTotal > vnl_math::eps ) gt = Gt[c] / GtTotal;
   			*(m_DataProbMapBuffer[c] + offset ) = gt;
   			postP[c] = ( gt * (Px[c] / Zx) );

			if (vnl_math_isnan(postP[c]) || vnl_math_isinf(postP[c])) {
				postP[c] = 0.0;
			}

			postPTotal += postP[c];
   		}

   		for ( size_t c = 0; c < nClasses; c++ ){
   			normP = ( postPTotal > 0.0 )?(postP[c]/postPTotal):uniformP;
			*(m_PostProbMapBuffer[c] + offset) = normP;
		}
	}
}

template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
void MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
SetSpatialEnergyMap( LabelType l, EnergyMapType* m ) {
	if ( m_SpatialEnergyMap.size() == 0 || m_SpatialEnergyMap.size() < m_NumberOfClasses ) {
		m_SpatialEnergyMap.resize( m_NumberOfClasses );
	}

	if ( l > 0 && l <= m_NumberOfClasses ) {
		float prob = 0.0;
		float energy = 0.0;
		size_t offset;
		size_t maxOffset = m->GetLargestPossibleRegion().GetNumberOfPixels();
		double corrFactor = ( 6 / -vcl_log( vnl_math::eps ) ) ;

		typename EnergyMapType::PixelType* eBuffer = m->GetBufferPointer();
	   	for ( offset = 0; offset < maxOffset; offset++) {
	   		prob = *( eBuffer + offset );
	   		energy = ( prob > vnl_math::eps )? ( -vcl_log( prob ) * this->m_ExternalLambda ) : MRF_INF ;
			*(m_SpatialEnergyMapBuffer[l-1] + offset) = energy ;
		}
	}
}

template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
void MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
SaveMRFEnergyMap( std::string name ) {
	typename EnergyMapType::Pointer map = EnergyMapType::New();
	map->SetRegions( m_Output->GetRequestedRegion() );
	map->CopyInformation( m_Output );
	map->Allocate();
	map->FillBuffer(-1.0);
	map->Update();

	EnergyValType* buffer = map->GetBufferPointer();

	size_t offset;
	size_t pos = 0;
	typename InputSampleType::ConstIterator end = m_Sample->End();
   	for (typename InputSampleType::ConstIterator it = m_Sample->Begin(); it != end; ++it, pos++) {
   		offset = it.GetInstanceIdentifier();
   		*( buffer + offset ) = GetCliqueEnergy( offset );
   	}

   	typedef itk::ImageFileWriter< EnergyMapType >  W;

   	typename W::Pointer w = W::New();
   	w->SetFileName( name );
   	w->SetInput( map );
   	w->Update();
}

template < class TInputImage, class TMaskImage, class TLabelType, class TEnergyVal >
void MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>::
SetModelParametersForLabel( LabelType l, ParametersType p ) {
	m_MembershipVector[l-1]->SetParameters( p );
	this->ModelModified();
}


#endif /* MRFENERGY_TXX_ */
