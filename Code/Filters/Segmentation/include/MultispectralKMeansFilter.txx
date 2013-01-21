// --------------------------------------------------------------------------------------
// File:    MultispectralKMeansFilter.txx
// Date:    Oct 4, 2011
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

#ifndef MULTISPECTRALKMEANSFILTER_TXX_
#define MULTISPECTRALKMEANSFILTER_TXX_

#include <itkTimeProbe.h>

#include "MultispectralKMeansFilter.h"

namespace itk {

/**
 *
 */
template <class TInputComponent, class TMaskImage>
MultispectralKMeansFilter< TInputComponent, TMaskImage>
::MultispectralKMeansFilter():
 	 m_NumberOfClasses(3), m_UserBucketSize(0), m_MaxIterations( 200 )
{
	// At least 1 inputs is necessary for a vector image.
	this->SetNumberOfRequiredInputs(0);
}

/**
 *
 */
template <class TInputComponent, class TMaskImage>
MultispectralKMeansFilter< TInputComponent, TMaskImage>
::~MultispectralKMeansFilter()
{
}


template <class TInputComponent, class TMaskImage>
void MultispectralKMeansFilter< TInputComponent, TMaskImage>
::Compute() {
	Superclass::GenerateData();
    itk::TimeProbe clock;

	m_NumberOfComponents = this->GetInputVectorImage()->GetNumberOfComponentsPerPixel();

	InputComponentPixelType minPixel = NumericTraits< InputComponentPixelType >::min();
	InputComponentPixelType maxPixel = NumericTraits< InputComponentPixelType >::max();

	MeasurementVectorType max( m_NumberOfComponents );
	max.Fill( minPixel );
	MeasurementVectorType min( m_NumberOfComponents );
	min.Fill( maxPixel );

	std::cout << "\t* Sample List Creation... ";
	std::cout.flush();
	clock.Start();

	typename InputVectorImageType::ConstPointer input = this->GetInputVectorImage();
	const typename InputVectorImageType::PixelContainer* buffer = input->GetPixelContainer();

	m_Sample = InputSampleType::New();
	m_Sample->SetImage( input );
	// Prepare mask for filters
	if( m_MaskImage.IsNotNull() ) {
		m_Sample->SetMask( m_MaskImage );
	}
	m_Sample->Initialize();


	size_t sampleSize = m_Sample->Size();

	if( m_NumberOfComponents!= m_Sample->GetMeasurementVectorSize() ) {
		itkExceptionMacro( "Input image number of components and samples mesasurement vector size don't match");
	}

	typename InputSampleType::ConstIterator end = m_Sample->End();

	for(typename InputSampleType::ConstIterator it = m_Sample->Begin(); it!=end; ++it) {
		for (size_t c = 0; c < m_NumberOfComponents; c++) {
			const double val = (*buffer)[it.GetInstanceIdentifier() + c];
			if ( val > max[c] ) max[c] = val;
			if ( val < min[c] ) min[c] = val;
		}
	}

	InputComponentPixelType absoluteMin = maxPixel;
	InputComponentPixelType absoluteMax = minPixel;

	for (size_t c = 0; c < m_NumberOfComponents; c++) {
		if ( max[c] > absoluteMax ) absoluteMax = max[c];
		if ( min[c] < absoluteMin ) absoluteMin = min[c];
	}
	clock.Stop();
	std::cout << "[DONE, mean=" << clock.GetMean() << " sec., total=" << clock.GetTotal() << " sec.]."<< std::endl;

    typename EstimatorType::ParametersType initialMeans( m_NumberOfClasses * m_NumberOfComponents );

    if( m_InitialParameters.size() == 0 ) {
    	initialMeans.Fill(0.0);

    	// Select initial means at random
    	for ( size_t cId = 0; cId < m_NumberOfClasses; cId++ ) {
    		MeasurementVectorType ci = m_Sample->GetRandomMeasurementVector();

    		for( size_t channel = 0; channel < m_NumberOfComponents; channel++ )
    			initialMeans[ cId*m_NumberOfComponents + channel ] = ci[channel];

    	}
    } else {
    	for ( size_t i = 0; i < m_NumberOfClasses; i++) {
    		for ( size_t j = 0; j < m_NumberOfComponents; j++) {
    			size_t idx=j*(m_NumberOfClasses)+i;
    			initialMeans[idx] = m_InitialParameters[i][j];
    		}
    	}

    }
    std::cout << "\t* Initial Means=" << initialMeans << std::endl;

	unsigned int bucketSize = (m_UserBucketSize)?m_UserBucketSize:(sampleSize * 0.001);
	std::cout << "\t* Tree generation (BucketSize=" << bucketSize << ")...";
	std::cout.flush();
	typename TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
    //treeGenerator->SetSample( m_Sample );
	treeGenerator->SetInputImage( input );
	treeGenerator->SetMaskImage( m_MaskImage );
    treeGenerator->SetBucketSize( bucketSize );


    clock.Start();
    treeGenerator->Update();
    clock.Stop();
    std::cout << "[DONE, mean=" << clock.GetMean() << ", total=" << clock.GetTotal() << "]."<< std::endl;

    m_Estimator = EstimatorType::New();
    m_Estimator->SetParameters( initialMeans );
    m_Estimator->SetKdTree( treeGenerator->GetOutput() );
    m_Estimator->SetMaximumIteration( m_MaxIterations );

    double centroidPositionThres = ( absoluteMax - absoluteMin ) * 1e-5;
    m_Estimator->SetCentroidPositionChangesThreshold( centroidPositionThres );

    std::cout << "\t* Estimator Optimization ( maxIterations= " << m_MaxIterations << ", centroidPositionChangesThres=" << centroidPositionThres << ") ...";
    std::cout.flush();

    clock.Start();
    m_Estimator->StartOptimization();
    clock.Stop();
    std::cout << "[DONE, mean=" << clock.GetMean() << ", total=" << clock.GetTotal() << "]."<< std::endl;

    m_FinalMeans = m_Estimator->GetParameters();

    std::cout << "\t* Final Means=" << m_FinalMeans << std::endl;

    for ( size_t l=0; l<m_NumberOfClasses; l++ ) {
    	OutputParametersType tissueMeansArray(m_NumberOfComponents);

    	for ( size_t c = 0; c< m_NumberOfComponents; c++ ) {
    		size_t idx = l*m_NumberOfComponents+c;
    		tissueMeansArray[c] = m_FinalMeans[idx];
    	}
    	m_OutputParameters.push_back( tissueMeansArray );
    }

    std::sort( m_OutputParameters.begin(), m_OutputParameters.end(), mySortMeans );
}


template <class TInputComponent, class TMaskImage>
void MultispectralKMeansFilter< TInputComponent, TMaskImage>
::GenerateData() {
	if( m_OutputParameters.size() == 0 ) this->Compute();

    typename DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();
    typename ImageClassifierType::Pointer classifier = ImageClassifierType::New();
    classifier->SetDecisionRule( (typename ImageClassifierType::DecisionRuleType*) decisionRule.GetPointer() );
    classifier->SetImage( this->GetInputVectorImage() );
    classifier->SetNumberOfClasses( m_NumberOfClasses );

    typename ImageClassLabelVectorObjectType::Pointer  classLabelsObject = ImageClassLabelVectorObjectType::New();
    classifier->SetClassLabels( classLabelsObject );

    ImageClassLabelVectorType &  classLabelsVector = classLabelsObject->Get();

    unsigned char firstLabel = ( m_MaskImage.IsNotNull() ) ? 1 : 0;
    for ( unsigned int i=0; i<m_NumberOfClasses; i++) {
    	unsigned char label = firstLabel + i;
    	classLabelsVector.push_back( label );
    }


    typename ImageMembershipFunctionVectorObjectType::Pointer membershipFunctionsObject = ImageMembershipFunctionVectorObjectType::New();
    classifier->SetMembershipFunctions( membershipFunctionsObject );

    ImageMembershipFunctionVectorType &  membershipFunctionsVector = membershipFunctionsObject->Get();

    typename MembershipFunctionType::CentroidType origin( m_NumberOfComponents );

    for ( unsigned int i = 0 ; i < m_NumberOfClasses ; i++ ) {
      MembershipFunctionPointer membershipFunction = MembershipFunctionType::New();
      for ( unsigned int j = 0 ; j < m_NumberOfComponents; j++ ) {
        origin[j] = m_OutputParameters[i][j];
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
    im_it.Begin();

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
void MultispectralKMeansFilter< TInputComponent, TMaskImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "ImageVector:" << std::endl;
  this->GetInputVectorImage()->Print( os, indent );
}

}

#endif /* MULTISPECTRALKMEANSFILTER_TXX_ */
