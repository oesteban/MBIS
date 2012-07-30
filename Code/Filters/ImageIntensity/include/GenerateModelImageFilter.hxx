// --------------------------------------------------------------------------------------
// File:          GenerateModelImageFilter.hxx
// Date:          May 21, 2012
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

#ifndef GENERATEMODELIMAGEFILTER_HXX_
#define GENERATEMODELIMAGEFILTER_HXX_

namespace mfbs {

/**
 *
 */
template <class TProbabilityInput, class TOutputPixelType >
GenerateModelImageFilter< TProbabilityInput, TOutputPixelType >
::GenerateModelImageFilter():
 m_UseBackgroundLabel(false)
{
  // At least 1 inputs is necessary for a vector image.
  this->SetNumberOfRequiredInputs( 0 );
}

template <class TProbabilityInput, class TOutputPixelType >
void GenerateModelImageFilter< TProbabilityInput, TOutputPixelType >
::GenerateData() {
	// Check that Paramters and Inputs match
	size_t numberOfLabels = m_Parameters.size();

	if( numberOfLabels < 1 || numberOfLabels!= m_InputVector.size() ) {
		itkExceptionMacro( << "Incorrect Number of Labels" );
	}

	size_t numberOfChannels =itk::NumericTraits< MeasurementVectorType >::GetLength( m_Parameters[0] );
	size_t numberOfPixels = m_InputVector[0]->GetLargestPossibleRegion().GetNumberOfPixels();

    // Generate empty data structure for output image
	MeasurementVectorType m( numberOfChannels );
	m.Fill( 0.0 );
    OutputImagePointer outputPtr = this->GetOutput(0);
    outputPtr->SetRegions( m_InputVector[0]->GetLargestPossibleRegion() );
    outputPtr->CopyInformation( m_InputVector[0] );
    outputPtr->SetNumberOfComponentsPerPixel( numberOfChannels );
    outputPtr->Allocate();
    outputPtr->FillBuffer( m );

    if( m_UseBackgroundLabel ) {
    	m_EstimatedBackground = ProbabilityImageType::New();
        m_EstimatedBackground->SetRegions( m_InputVector[0]->GetLargestPossibleRegion() );
        m_EstimatedBackground->CopyInformation( m_InputVector[0] );
        m_EstimatedBackground->Allocate();
        m_EstimatedBackground->FillBuffer( 0.0 );
    }


    //MeasurementVectorType* outBuffer = outputPtr->GetBufferPointer();

    std::vector< const ProbabilityPixelType* > probBuffer;
    for( size_t i = 0; i<numberOfLabels; i++)
    	probBuffer.push_back( m_InputVector[i]->GetBufferPointer() );


    ProbabilityPixelType p = 0.0;
    ProbabilityPixelType totalP;
    ProbabilityPixelType maxP;
    size_t mapLabel;

	for( size_t pixId = 0; pixId < numberOfPixels; pixId++ ) {
		maxP = 0.0;
		mapLabel = 0;
		m.Fill(0.0);
		totalP = 0.0;

		for ( size_t l = 0; l < numberOfLabels; l++  ){
			p = *(probBuffer[l] + pixId);
			if( p > maxP ) {
				maxP = p;
				mapLabel = l;
			}
			totalP+=p;

			for( size_t c = 0; c < numberOfChannels; c++ ) {
				m[c] +=  p * m_Parameters[l][c];
			}
		}

		if( totalP > 0.0 ) {
			outputPtr->SetPixel( outputPtr->ComputeIndex(pixId), m );
		//*(outBuffer + pixId*numberOfChannels ) = val;
		}

		if( m_UseBackgroundLabel && mapLabel!=0 ) {
			m_EstimatedBackground->SetPixel( outputPtr->ComputeIndex(pixId), 1.0 );
		}
	}



}

}


#endif /* GENERATEMODELIMAGEFILTER_HXX_ */
