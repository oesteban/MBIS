// --------------------------------------------------------------------------------------
// File:          MultispectralGaussianEstimator.txx
// Date:          Nov 14, 2011
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

#ifndef MULTISPECTRALGAUSSIANESTIMATOR_TXX_
#define MULTISPECTRALGAUSSIANESTIMATOR_TXX_

#include "MultispectralGaussianEstimator.h"

#include <itkLabelMapMaskImageFilter.h>
#include <itkLabelImageToLabelMapFilter.h>

namespace itk {

/**
 *
 */
template <class TInputComponent, class TProbabilityPixelType>
MultispectralGaussianEstimator< TInputComponent, TProbabilityPixelType >
::MultispectralGaussianEstimator():
 m_MaximumIteration(3), m_NumberOfClasses(3), m_UseMaskedMode(true), m_UsePriorProbabilityImages(false)
{
	// At least 1 inputs is necessary for a vector image.
	m_Priors.resize(m_NumberOfClasses);
	this->SetNumberOfRequiredInputs(0);
}

/**
 *
 */
template <class TInputComponent, class TProbabilityPixelType>
MultispectralGaussianEstimator< TInputComponent, TProbabilityPixelType >
::~MultispectralGaussianEstimator()
{
}


template <class TInputComponent, class TProbabilityPixelType>
void MultispectralGaussianEstimator< TInputComponent, TProbabilityPixelType >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "ImageVector:" << std::endl;
  m_Input->Print( os, indent );
}

template <class TInputComponent, class TProbabilityPixelType>
void MultispectralGaussianEstimator< TInputComponent, TProbabilityPixelType >
::SetPriorProbabilityImage( unsigned int classId,
		    const MultispectralGaussianEstimator< TInputComponent, TProbabilityPixelType >::ProbabilityImageType * image )
{
	if ( classId >= m_NumberOfClasses )
		itkExceptionMacro("Trying to set a prior probability image for a non-existant class (id " << classId << "> lastClassId " << (m_NumberOfClasses-1) << ").");
	m_Priors[classId]= image;
}


template <class TInputComponent, class TProbabilityPixelType>
void MultispectralGaussianEstimator< TInputComponent, TProbabilityPixelType >
::SetPriors( const MultispectralGaussianEstimator< TInputComponent, TProbabilityPixelType >::ProbabilityImagesVector& priors ) {
	if( priors.size() != m_NumberOfClasses )
		itkExceptionMacro("Number of priors and classes do not match.");

	typename ProbabilityImagesVector::const_iterator it = priors.beguin();
	typename ProbabilityImagesVector::const_iterator end = priors.end();

	unsigned int classId = 0;
	while( it!= end ) {
		m_Priors[classId] = *it;
		classId++;
		it++;
	}
}

// Compute Weighted Mean & Cov.
// http://www.itk.org/Doxygen/html/classitk_1_1Statistics_1_1WeightedMeanSampleFilter.html
// http://www.itk.org/Doxygen/html/classitk_1_1Statistics_1_1WeightedCovarianceSampleFilter.html
template <class TInputComponent, class TProbabilityPixelType>
void MultispectralGaussianEstimator< TInputComponent, TProbabilityPixelType >
::GenerateData() {
	m_Input = this->GetInputVectorImage();
    m_NumberOfComponents = m_Input->GetNumberOfComponentsPerPixel();
	m_PriorIndexOffset = (unsigned int) m_MaskImage.IsNotNull();

	m_Sample = InputSampleType::New();
	m_Sample->SetImage( m_Input );

	// Check the agreement between priors and classes
	if( m_Priors.size() != m_NumberOfClasses ) {
		itkExceptionMacro( "Insufficient prior specified: " << m_Priors.size() << " when " << (int) m_NumberOfClasses << " are required.");
	}

	// Prepare mask for filters
	ProbabilitySamplePointer maskSample = ProbabilitySampleType::New();
	maskSample->SetImage( m_MaskImage );
	maskSample->Update();

	// Weighted Mean & Covariance estimators
	typename CovarianceEstimatorType::Pointer covEst = CovarianceEstimatorType::New();
	covEst->SetInput(m_Sample);
	const typename CovarianceEstimatorType::MatrixDecoratedType * cov_out = covEst->GetCovarianceMatrixOutput();
	const typename CovarianceEstimatorType::MeasurementVectorDecoratedType * mean_out = covEst->GetMeanOutput();

	// Compute Mean & Covariance matrix for each class
	for ( unsigned int i = 0; i<m_NumberOfClasses; i++) {
		// Weights array generation for this class
		ProbabilitySamplePointer wSample = ProbabilitySampleType::New();
		wSample->SetImage( m_Priors[i] );

		WeightArrayType weights( wSample->Size() );
		weights.Fill(1.0);

		typename ProbabilitySampleType::ConstIterator end = wSample->End();
		for(typename ProbabilitySampleType::ConstIterator it = wSample->Begin(); it!=end; ++it ) {
			unsigned int iid = it.GetInstanceIdentifier();
			weights[iid] = it.GetMeasurementVector()[0];
			(m_CVector[i])->SetWeight( iid, weights[iid] );
		}

		covEst->SetWeights( weights );

		try {
			covEst->Update();
		} catch( itk::ExceptionObject & ex ) {
			std::cerr << "Exception caught" << ex << std::endl;
		}

		ParametersType params( (m_NumberOfComponents +1 ) * m_NumberOfComponents);

		for ( unsigned int k = 0; k < mean_out->Get().Size(); k++)
			params[k] = mean_out->Get()[k];

		for ( unsigned int j = 0; j < cov_out->Get().Rows(); j++)
			for (unsigned int k = 0; k< cov_out->Get().Cols(); k++)
				params[(j+1)*m_NumberOfComponents + k] = cov_out->Get()(k,j);

		std::cout << "\tClass [" << i << "]: " << params << std::endl;

		m_StatParameters.push_back( params );

	}

}

template <class TInputComponent, class TProbabilityPixelType>
void MultispectralGaussianEstimator< TInputComponent, TProbabilityPixelType >
::SetMeans( ParametersVectorType& means ) {

	ProbabilityImagePointer unoImage = ProbabilityImageType::New();
	unoImage->SetRegions( this->GetInputVectorImage()->GetLargestPossibleRegion() );
	unoImage->CopyInformation( this->GetInputVectorImage() );
	unoImage->Allocate();
	unoImage->Fill( 1.0 );
	unoImage->Update();

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

    for ( size_t i = 0 ; i < m_NumberOfClasses ; i++ ) {
      MembershipFunctionPointer membershipFunction = MembershipFunctionType::New();
      for ( size_t j = 0 ; j < m_NumberOfComponents; j++ ) {
        origin[j] = means[i][j];
      }
      membershipFunction->SetCentroid( origin );
      membershipFunctionsVector.push_back( membershipFunction.GetPointer() );
    }

    classifier->Update();

    OutputImagePointer out = classifier->GetOutput();


	typedef itk::LabelImageToLabelMapFilter< OutputImageType > ClassifiedToLabelMapFilter;
	typedef typename ClassifiedToLabelMapFilter::OutputImageType LabelMap;
	typedef itk::LabelMapMaskImageFilter< LabelMap, ProbabilityImageType > LabelMapMaskFilter;

    for ( typename OutputImageType::PixelType i = 0 ; i < m_NumberOfClasses ; i++ ) {
		typename ClassifiedToLabelMapFilter::Pointer lfilter = ClassifiedToLabelMapFilter::New();
		lfilter->SetInput( out->GetOutput() );
		lfilter->Update();

		typename LabelMapMaskFilter::Pointer lmask = LabelMapMaskFilter::New();
		lmask->SetInput1( lfilter->GetOutput() );
		lmask->SetInput2( unoImage );
		lmask->SetLabel( i + firstLabel );
		lmask->Update();

		this->SetPriorProbabilityImage( i, lmask->GetOutput());
    }
}

}


#endif /* MULTISPECTRALGAUSSIANESTIMATOR_TXX_ */
