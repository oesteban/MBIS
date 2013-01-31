// --------------------------------------------------------------------------------------
// File:          NormalizeToModelImageFilter.hxx
// Date:          Jan 30, 2013
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2013, code@oscaresteban.es (Oscar Esteban)
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

#ifndef NORMALIZETOMODELIMAGEFILTER_HXX_
#define NORMALIZETOMODELIMAGEFILTER_HXX_

#include "NormalizeToModelImageFilter.h"

namespace mfbs {

/**
 *
 */
template<class TInputVectorImage >
NormalizeToModelImageFilter< TInputVectorImage >
::NormalizeToModelImageFilter(): m_Lambda(1.0) {
  // At least 1 inputs is necessary for a vector image.
  this->SetNumberOfRequiredInputs( 2 );
}

template<class TInputVectorImage >
void NormalizeToModelImageFilter< TInputVectorImage >
::GenerateData() {
	size_t numChannels = this->GetInput(0)->GetNumberOfComponentsPerPixel();
	size_t numPixels = this->GetInput(0)->GetLargestPossibleRegion().GetNumberOfPixels();
	bool isMasked = m_MaskImage.IsNotNull();

	const float* maskBuffer;

	if ( isMasked ) {
		maskBuffer = m_MaskImage->GetBufferPointer();
	}

	// Generate empty data structure for output image
	InputVectorImagePointer outputPtr = this->GetOutput(0);
	outputPtr->SetRegions( this->GetInput(0)->GetRequestedRegion() );
	outputPtr->CopyInformation( this->GetInput(0) );
	outputPtr->SetNumberOfComponentsPerPixel( numChannels );
	outputPtr->Allocate();
	MeasurementVectorType m ( numChannels );
	m.Fill(0.0);
	outputPtr->FillBuffer( m );

	typename OutputImageType::InternalPixelType* outBuffer = outputPtr->GetBufferPointer();

	typename SelectionFilter::Pointer sel0 = SelectionFilter::New();
	sel0->SetInput( this->GetInput(0) );

	typename SelectionFilter::Pointer sel1 = SelectionFilter::New();
	sel1->SetInput( this->GetInput(1) );

	std::vector< MeasurementType > sample0;
	std::vector< MeasurementType > sample1;

	for ( size_t i = 0; i < numChannels; i++ ){
		// Get component i
		sel0->SetIndex(i);
		sel1->SetIndex(i);
		sel0->Update();
		sel1->Update();

		// Compute max, min
		const MeasurementType* buf0 = sel0->GetOutput()->GetBufferPointer();
		const MeasurementType* buf1 = sel1->GetOutput()->GetBufferPointer();

		for( size_t p = 0; p<numPixels; p++ ) {
			if( !isMasked || *( maskBuffer+p) > 0 ){
				sample0.push_back( *( buf0 + p ) );
				sample1.push_back( *( buf1 + p ) );
			}
		}

		std::sort( sample0.begin(), sample0.end() );
		std::sort( sample1.begin(), sample1.end() );
		size_t sampleSize = sample0.size();
		MeasurementType maximum[2],minimum[2];
		maximum[0] = sample0[ (size_t) (sampleSize-1)*0.95 ];
        maximum[1] = sample1[ (size_t) (sampleSize-1)*0.95 ];
		minimum[0] = sample0[ (size_t) (sampleSize-1)*0.05 ];
		minimum[1] = sample1[ (size_t) (sampleSize-1)*0.05 ];
		MeasurementType mean2 = sample1[ (size_t) (sampleSize-1)*0.5 ];

		MeasurementType A = (maximum[0]-minimum[0]) / (maximum[1]-minimum[1]);
		//MeasurementType B = minimum[0] - minimum[1] * A;
		MeasurementType B = - mean2*A;

		// Normalize & set onto bias field vector
		MeasurementVectorType m;
		for ( size_t p = 0; p< numPixels; p++) {
			typename OutputImageType::IndexType idx = outputPtr->ComputeIndex( p );
			if( !isMasked || *( maskBuffer+p) > 0 ){
				m = outputPtr->GetPixel( idx );
				m[i] = m_Lambda * (*( buf1 + p ) * A + B );
				outputPtr->SetPixel( idx, m  );
			}
		}
	}



}

}

#endif /* NORMALIZETOMODELIMAGEFILTER_HXX_ */
