// --------------------------------------------------------------------------------------
// File:          BSplineRegressionImageFilter.hxx
// Date:          May 18, 2012
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

#ifndef BSPLINEREGRESSIONIMAGEFILTER_HXX_
#define BSPLINEREGRESSIONIMAGEFILTER_HXX_

//#include "Common/include/itkImportImageFilter.h"
#include <itkBSplineControlPointImageFilter.h>
#include "BSplineRegressionImageFilter.h"
#include <itkMedianImageFilter.h>

namespace mfbs {


template< class TInputVectorImage >
BSplineRegressionImageFilter< TInputVectorImage >
::BSplineRegressionImageFilter() :
   m_SplineOrder( 3 )
{
	this->SetNumberOfRequiredInputs( 1 );

	this->m_NumberOfControlPoints.Fill( 4 );
}

template< class TInputVectorImage >
void
BSplineRegressionImageFilter< TInputVectorImage >
::GenerateData() {
	size_t numChannels = this->GetInput(0)->GetNumberOfComponentsPerPixel();

	bool isMasked = m_MaskImage.IsNotNull();

	// Generate empty data structure for output image
	InputVectorImagePointer outputPtr = this->GetOutput(0);
	outputPtr->SetRegions( this->GetInput(0)->GetRequestedRegion() );
	outputPtr->CopyInformation( this->GetInput(0) );
	outputPtr->SetNumberOfComponentsPerPixel( numChannels );
	outputPtr->Allocate();
	MeasurementVectorType m ( numChannels );
	m.Fill(0.0);
	outputPtr->FillBuffer( m );


	typename BSplineFilterType::ArrayType numberOfFittingLevels;
	numberOfFittingLevels.Fill( 1 );

	typename ScalarImageType::PointType parametricOrigin =  this->GetInput(0)->GetOrigin();
	for( unsigned int d = 0; d < ImageDimension; d++ ) {
	  parametricOrigin[d] += (
			  this->GetInput(0)->GetSpacing()[d] *
			  this->GetInput(0)->GetLargestPossibleRegion().GetIndex()[d] );
	}


	// Temporarily set the direction cosine to identity since the B-spline
	// approximation algorithm works in parametric space and not physical
	// space.
	typename ScalarImageType::DirectionType identity;
	identity.SetIdentity();

	const typename ScalarImageType::RegionType & bufferedRegion = this->GetInput(0)->GetBufferedRegion();
	const SizeValueType numberOfPixels = bufferedRegion.GetNumberOfPixels();
	const bool filterHandlesMemory = false;

	// Initialize a copy of input
	typedef itk::ImageDuplicator< InputVectorImageType > Duplicator;
	typename Duplicator::Pointer d = Duplicator::New();
	d->SetInputImage( this->GetInput(0) );
	d->Update();
	InputVectorImagePointer input = d->GetOutput();
	input->SetDirection( identity );

	for ( size_t channel = 0; channel < numChannels; channel++ ) {
		// Split input into its channels
		//typedef VectorIndexSelectionCastImageFilter< typename ImporterType::OutputImageType, ChannelImageType> NewSelectorType;
		typename InputSelectorType::Pointer channelSelector = InputSelectorType::New();
		channelSelector->SetInput( input );
		channelSelector->SetIndex( channel );
		channelSelector->Update();
		channelSelector->GetOutput()->SetRegions( input->GetRequestedRegion() );
		typename ChannelImageType::Pointer channelEstimate = channelSelector->GetOutput();

		PointSetPointer fieldPoints = PointSetType::New();
		fieldPoints->Initialize();

		typename BSplineFilterType::WeightsContainerType::Pointer weights =  BSplineFilterType::WeightsContainerType::New();
		weights->Initialize();

		ImageRegionConstIteratorWithIndex<ChannelImageType> It( channelEstimate, channelEstimate->GetRequestedRegion() );

		unsigned int index = 0;

		ScalarType scalar;
		PointType point;
		for ( It.GoToBegin(); !It.IsAtEnd(); ++It ) {
			channelEstimate->TransformIndexToPhysicalPoint( It.GetIndex(), point );
			if ( !isMasked || m_MaskImage->GetPixel( It.GetIndex() ) > 0 ) {
				scalar[0] = It.Get();
				weights->InsertElement( index, 1.0 );
				fieldPoints->SetPointData( index, scalar );
				fieldPoints->SetPoint( index, point );

				index++;
			}

			++It;
			if ( It.IsAtEnd() ) break;
			++It;
			if ( It.IsAtEnd() ) break;
			++It;
			if ( It.IsAtEnd() ) break;
		}

		typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
		bspliner->SetOrigin( parametricOrigin );
		bspliner->SetSpacing( channelEstimate->GetSpacing() );
		bspliner->SetSize( channelEstimate->GetLargestPossibleRegion().GetSize() );
		bspliner->SetDirection( channelEstimate->GetDirection() );
		bspliner->SetGenerateOutputImage( false );
		bspliner->SetNumberOfLevels( numberOfFittingLevels );
		bspliner->SetSplineOrder( this->m_SplineOrder );
		bspliner->SetNumberOfControlPoints( m_NumberOfControlPoints );
		bspliner->SetInput( fieldPoints );
		bspliner->SetPointWeights( weights );
		bspliner->Update();

		// Add the bias field control points to the current estimate.
		typename BiasFieldControlPointLatticeType::Pointer controlPointLattice = bspliner->GetPhiLattice();

		typedef BSplineControlPointImageFilter
		<BiasFieldControlPointLatticeType, ScalarImageType> BSplineReconstructerType;
		typename BSplineReconstructerType::Pointer reconstructer =
		  BSplineReconstructerType::New();
		reconstructer->SetInput( controlPointLattice );
		reconstructer->SetOrigin( channelEstimate->GetOrigin() );
		reconstructer->SetSpacing( channelEstimate->GetSpacing() );
		reconstructer->SetDirection( channelEstimate->GetDirection() );
		reconstructer->SetSize( channelEstimate->GetLargestPossibleRegion().GetSize() );
		reconstructer->Update();

		typedef VectorIndexSelectionCastImageFilter<ScalarImageType, ChannelImageType> SelectorType;
		typename SelectorType::Pointer selector = SelectorType::New();
		selector->SetInput( reconstructer->GetOutput() );
		selector->SetIndex( 0 );
		selector->Update();
		selector->GetOutput()->SetRegions( channelEstimate->GetRequestedRegion() );

		typename ChannelImageType::Pointer smoothField = selector->GetOutput();
		smoothField->Update();
		smoothField->DisconnectPipeline();
		smoothField->SetRegions( channelEstimate->GetRequestedRegion() );

		// Combine to bias field vector
		for ( size_t i = 0; i< numberOfPixels; i++) {
			typename InputVectorImageType::IndexType idx = smoothField->ComputeIndex( i );
			if ( !isMasked || m_MaskImage->GetPixel( idx ) > 0 ) {
				m = outputPtr->GetPixel( idx );
				m[channel] = smoothField->GetPixel( idx );
				outputPtr->SetPixel( idx, m );
			}
		}


	}

}

template< class TInputVectorImage >
void
BSplineRegressionImageFilter< TInputVectorImage >
::PrintSelf( std::ostream &os, itk::Indent indent ) const {


}

}


#endif /* BSPLINEREGRESSIONIMAGEFILTER_HXX_ */
