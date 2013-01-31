// --------------------------------------------------------------------------------------
// File:          MultispectralMLClassifierFilter.hxx
// Date:          Jan 24, 2013
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

#ifndef MULTISPECTRALMLCLASSIFIERFILTER_HXX_
#define MULTISPECTRALMLCLASSIFIERFILTER_HXX_

#include "MultispectralMLClassifierFilter.h"

namespace itk {

/**
 *
 */
template <class TInputComponent, class TMaskImage>
MultispectralMLClassifierFilter< TInputComponent, TMaskImage>
::MultispectralMLClassifierFilter() {
	// At least 1 inputs is necessary for a vector image.
	this->SetNumberOfRequiredInputs(0);
}


template <class TInputComponent, class TMaskImage>
void MultispectralMLClassifierFilter< TInputComponent, TMaskImage>
::GenerateData() {
	Superclass::GenerateData();


	if( m_Parameters.size() == 0 ) {
		itkExceptionMacro( "Parameters are not set" );
	}

	std::sort( m_Parameters.begin(), m_Parameters.end(), mySortMeans );
	size_t numberOfClasses = m_Parameters.size();
	size_t numberOfComponents = this->GetInputVectorImage()->GetNumberOfComponentsPerPixel();

    typename DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();
    typename ImageClassifierType::Pointer classifier = ImageClassifierType::New();
    classifier->SetDecisionRule( (typename ImageClassifierType::DecisionRuleType*) decisionRule.GetPointer() );
    classifier->SetImage( this->GetInputVectorImage() );
    classifier->SetNumberOfClasses( numberOfClasses );

    typename ImageClassLabelVectorObjectType::Pointer  classLabelsObject = ImageClassLabelVectorObjectType::New();
    classifier->SetClassLabels( classLabelsObject );

    ImageClassLabelVectorType &  classLabelsVector = classLabelsObject->Get();

    unsigned char firstLabel = ( m_MaskImage.IsNotNull() ) ? 1 : 0;
    for ( unsigned int i=0; i<numberOfClasses; i++) {
    	unsigned char label = firstLabel + i;
    	classLabelsVector.push_back( label );
    }


    typename ImageMembershipFunctionVectorObjectType::Pointer membershipFunctionsObject = ImageMembershipFunctionVectorObjectType::New();
    classifier->SetMembershipFunctions( membershipFunctionsObject );

    ImageMembershipFunctionVectorType &  membershipFunctionsVector = membershipFunctionsObject->Get();

    typename MembershipFunctionType::CentroidType origin( numberOfComponents );

    for ( unsigned int i = 0 ; i < numberOfClasses ; i++ ) {
      MembershipFunctionPointer membershipFunction = MembershipFunctionType::New();
      for ( unsigned int j = 0 ; j < numberOfComponents; j++ ) {
        origin[j] = m_Parameters[i][j];
      }
      membershipFunction->SetCentroid( origin );
      membershipFunctionsVector.push_back( membershipFunction.GetPointer() );
    }

    classifier->Update();

    OutputImagePointer out = classifier->GetOutput();

    // Generate empty data structure for output image
    OutputImagePointer outputPtr = this->GetOutput(0);
    outputPtr->SetRegions( this->GetInputVectorImage()->GetLargestPossibleRegion() );
    outputPtr->SetSpacing( this->GetInputVectorImage()->GetSpacing() );
    outputPtr->SetDirection( this->GetInputVectorImage()->GetDirection() );
    outputPtr->SetOrigin( this->GetInputVectorImage()->GetOrigin() );
    outputPtr->Allocate();
    outputPtr->FillBuffer( 0.0 );

    // Iterator to copy to output
    OutputImageIterator im_it( outputPtr, outputPtr->GetLargestPossibleRegion());
    im_it.GoToBegin();

    // Generate (copy data to output). Do not copy pixels outside masks or outside thresholds window (min1,max1)
    while ( !im_it.IsAtEnd() ) {
    	if (  m_MaskImage.IsNull() ||
    	     (m_MaskImage.IsNotNull() && m_MaskImage->GetPixel(im_it.GetIndex() )) ) {
    		im_it.Set( out->GetPixel( im_it.GetIndex() ) );
    	}
    	++im_it;
    }
}

template <class TInputComponent, class TMaskImage>
void MultispectralMLClassifierFilter< TInputComponent, TMaskImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "ImageVector:" << std::endl;
  this->GetInputVectorImage()->Print( os, indent );
}

}


#endif /* MULTISPECTRALMLCLASSIFIERFILTER_HXX_ */
