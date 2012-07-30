// --------------------------------------------------------------------------------------
// File:          ImageBiasedGMModelEMEstimator.hxx
// Date:          May 14, 2012
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

#ifndef IMAGEBIASEDGMMODELEMESTIMATOR_HXX_
#define IMAGEBIASEDGMMODELEMESTIMATOR_HXX_

#include <numeric>
#include <itkNumericTraits.h>
#include <itkImageDuplicator.h>
#include "ImageBiasedGMModelEMEstimator.h"


namespace mfbs {
namespace Statistics {
template <class TInputVectorImage, class TProbabilityPixelType>
ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::ImageBiasedGMModelEMEstimator() {
	m_TerminationCode = NOT_CONVERGED;

	m_MembershipFunctionsObject = MembershipFunctionVectorObjectType::New();
	m_MembershipFunctionsWeightArrayObject = MembershipFunctionsWeightsArrayObjectType::New();
	m_Sample = 0;
	m_MaxIteration = 100;
	m_UseBiasCorrection = true;
}

template <class TInputVectorImage, class TProbabilityPixelType>
void ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::PrintSelf(std::ostream & os, Indent indent) const {
	Superclass::PrintSelf(os, indent);
	os << indent << "Maximum Iteration: " << this->GetMaximumIteration() << std::endl;
	os << indent << "Sample: " << m_Sample << std::endl;
	os << indent << "Number Of Components: " << this->GetNumberOfComponents() << std::endl;
	for (unsigned int i = 0; i < this->GetNumberOfComponents(); i++) {
		os << indent << "Component Membership Function[" << i << "]: " << this->GetComponentMembershipFunction(i) << std::endl;
	}
	os << indent << "Termination Code: " << this->GetTerminationCode() << std::endl;
	os << indent << "Initial Proportions: " << this->GetInitialProportions() << std::endl;
	os << indent << "Proportions: " << this->GetProportions() << std::endl;
	os << indent << "Calculated Expectation: " << this->CalculateExpectation() << std::endl;
}

template <class TInputVectorImage, class TProbabilityPixelType>
void ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::SetMaximumIteration(int numberOfIterations) {
	m_MaxIteration = numberOfIterations;
}

template <class TInputVectorImage, class TProbabilityPixelType>
int ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::GetMaximumIteration() const {
	return m_MaxIteration;
}

template <class TInputVectorImage, class TProbabilityPixelType>
void ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::SetInitialProportions(ProportionVectorType & proportions) {
	m_InitialProportions = proportions;
}

template <class TInputVectorImage, class TProbabilityPixelType>
const typename ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::ProportionVectorType &
ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::GetInitialProportions() const {
	return m_InitialProportions;
}

template <class TInputVectorImage, class TProbabilityPixelType>
const typename ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::ProportionVectorType &
ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::GetProportions() const {
	return m_Proportions;
}

template <class TInputVectorImage, class TProbabilityPixelType>
int ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::AddComponent(ComponentType *component) {
	m_ComponentVector.push_back(component);
	return static_cast<int>(m_ComponentVector.size());
}

template <class TInputVectorImage, class TProbabilityPixelType>
unsigned int ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::GetNumberOfComponents() const {
	return m_ComponentVector.size();
}

template <class TInputVectorImage, class TProbabilityPixelType>
typename ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::TERMINATION_CODE
ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::GetTerminationCode() const {
	return m_TerminationCode;
}

template <class TInputVectorImage, class TProbabilityPixelType>
typename ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::ComponentMembershipFunctionType *
ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::GetComponentMembershipFunction(int componentIndex) const {
	return (m_ComponentVector[componentIndex])->GetMembershipFunction();
}

template <class TInputVectorImage, class TProbabilityPixelType>
bool ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::CalculateDensities() {
	// Formula (25) in Bach2005. A posteriori probabilities
	bool componentModified = false;

	for (size_t i = 0; i < m_ComponentVector.size(); i++) {
		if ((m_ComponentVector[i])->AreParametersModified()) {
			componentModified = true;
			break;
		}
	}

	if (!componentModified) {
		return false;
	}

	double temp;
	size_t numberOfComponents = m_ComponentVector.size();
	std::vector<double> tempWeights(numberOfComponents, 0.);
	double uniformWeight = 1.0 / numberOfComponents;

	typename InputSampleType::ConstIterator last = m_Sample->End();

	typedef typename InputSampleType::AbsoluteFrequencyType FrequencyType;
	FrequencyType frequency;
	FrequencyType zeroFrequency = NumericTraits < FrequencyType > ::Zero;
	typename InputSampleType::MeasurementVectorType yVector;
	double density;
	double densitySum;
	double totalDensitySum;
	double minDouble = NumericTraits<double>::epsilon();

	SizeValueType measurementVectorIndex = 0;

	for (typename InputSampleType::ConstIterator iter = m_Sample->Begin(); iter != last; ++iter, ++measurementVectorIndex) {
		yVector = iter.GetMeasurementVector();
		densitySum = 0.0;

		// For every x belonging Lpm
		for (size_t x = 0; x < numberOfComponents; ++x) {
			double Px = m_Proportions[x];
			double Py_x = m_ComponentVector[x]->Evaluate(yVector);

			// This is P(x) * P(y=yVector|x, theta )
			tempWeights[x] = Px * Py_x;
		}

		// This is the denominator of (25) in Bach2005
		for (size_t x = 0; x < numberOfComponents; ++x)
			densitySum += tempWeights[x];

		for (size_t x = 0; x < numberOfComponents; ++x) {
			temp = tempWeights[x];

			// Final division to get P(k)(x|y_i,theta)
			if (densitySum > minDouble) // just to make sure temp does not blow up!
				temp /= densitySum;
			else {
				temp = uniformWeight;
			}

			// Save weight
			// m_Ptly[measurementVectorIndex][x] = temp;
			*(m_Posteriors[x]->GetBufferPointer() + iter.GetInstanceIdentifier()) = temp;
		}
	}


	return true;
}

template <class TInputVectorImage, class TProbabilityPixelType>
double ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::CalculateExpectation() const {
	double sum = 0.0;

	if (m_Sample) {
		unsigned int measurementVectorIndex;
		SizeValueType size = m_Sample->Size();
		double logProportion;
		double temp;
		for (size_t componentIndex = 0; componentIndex < m_ComponentVector.size(); ++componentIndex) {
			temp = m_Proportions[componentIndex];

			// if temp is below the smallest positive double number
			// the log may blow up
			if (temp > NumericTraits<double>::epsilon()) {
				logProportion = vcl_log(temp);
			} else {
				logProportion = NumericTraits<double>::NonpositiveMin();
			}

			for (measurementVectorIndex = 0; measurementVectorIndex < size; measurementVectorIndex++) {
				temp = m_ComponentVector[componentIndex]->GetWeight(measurementVectorIndex);

				if (temp > NumericTraits<double>::epsilon()) {
					sum += temp * (logProportion + vcl_log(temp));
				} else {
					// let's throw an exception
                    itkExceptionMacro( << "temp is null" );
				}
			}
		}
	}
	return sum;
}

template <class TInputVectorImage, class TProbabilityPixelType>
bool ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::UpdateComponentParameters() {
	bool updated = false;
	size_t numberOfComponents = m_ComponentVector.size();
	ComponentType *component;

	typename InputSampleType::ConstIterator iter = m_Sample->Begin();
	typename InputSampleType::ConstIterator last = m_Sample->End();

	size_t i = 0;
	double totalWeight;
	double temp;
	double minDouble = NumericTraits<double>::epsilon();

	while (iter != last) {
		totalWeight = 0.;
		size_t idx = iter.GetInstanceIdentifier();

		for (size_t x = 0; x < numberOfComponents; ++x) {
			totalWeight += *(m_Posteriors[x]->GetBufferPointer() + idx );
			//totalWeight += m_Ptly[i][x];
		}

		for (size_t x = 0; x < numberOfComponents; ++x) {
			temp = 0.;
			if (totalWeight > minDouble) // just to make sure temp does not blow up!
				temp = (*(m_Posteriors[x]->GetBufferPointer() + idx )) / totalWeight;
				//temp = m_Ptly[i][x] / totalWeight;

			m_ComponentVector[x]->SetWeight(i, temp);
		}
		++iter;
		i++;
	}

	for (size_t x = 0; x < numberOfComponents; ++x) {
		component = m_ComponentVector[x];
		component->Update();
		updated = updated || component->AreParametersModified();
	}

	return updated;
}

template <class TInputVectorImage, class TProbabilityPixelType>
bool ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::UpdateProportions() {
	// TODO Compute and cache the histogram so this is computed like SUM( h(bin) * GetWeight(j) );
	size_t numberOfComponents = m_ComponentVector.size();
	size_t sampleSize = m_Sample->Size();
	double totalFrequency = static_cast<double>(m_Sample->GetTotalFrequency());
	bool updated = false;

	if (totalFrequency < NumericTraits<double>::epsilon()) {
		return updated;
	}

	std::vector<double> tempSum(numberOfComponents,0.0);
	typename InputSampleType::ConstIterator last = m_Sample->End();

	for (typename InputSampleType::ConstIterator iter = m_Sample->Begin(); iter != last; ++iter) {
		size_t idx = iter.GetInstanceIdentifier();

		for (size_t x = 0; x < numberOfComponents; ++x) {
			double freq= (*(m_Posteriors[x]->GetBufferPointer() + idx )) * iter.GetFrequency();
			tempSum[x] += freq;
		}
	}

	double totalSum = std::accumulate(tempSum.begin(),tempSum.end(),0);


	for (size_t i = 0; i < numberOfComponents; ++i) {
		tempSum[i] /= totalSum;

		if (tempSum[i] != m_Proportions[i]) {
			m_Proportions[i] = tempSum[i];
			updated = true;
		}
	}

	return updated;
}

template <class TInputVectorImage, class TProbabilityPixelType>
bool ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::UpdateBiasFieldEstimate() {
	size_t numberOfComponents = m_ComponentVector.size();
	//size_t biasOffset = ( m_MaskImage.IsNull() )?1:0;
	size_t numberOfChannels = m_Input->GetNumberOfComponentsPerPixel();
	//size_t numberOfBiasedLabels = numberOfComponents - biasOffset;

	// Get means
	std::vector< MeasurementVectorType > means;
	means.resize( numberOfComponents );

	for ( size_t i = 0; i < numberOfComponents; i++  ){
		itk::NumericTraits< MeasurementVectorType >::SetLength( means[i], numberOfChannels );

		typename ComponentType::ParametersType params = m_ComponentVector[i]->GetFullParameters();
		for( size_t j = 0; j < numberOfChannels; j++ ) {
			means[i][j] = params[j];
		}
	}

	// Generate current modelled image
	typename GenerateModelFilter::Pointer modelFilter = GenerateModelFilter::New();
	modelFilter->SetParameters( means );
	for(size_t i = 0; i< numberOfComponents; i++)
		modelFilter->AddInput( m_Posteriors[i] );

	modelFilter->SetUseBackgroundLabel( m_MaskImage.IsNull() );
	modelFilter->Update();

	// Log both images
	typename LogFilter::Pointer logModel = LogFilter::New();
	logModel->SetInput( modelFilter->GetOutput() );
	logModel->Update();

	typename LogFilter::Pointer logInput = LogFilter::New();
	logInput->SetInput( m_CorrectedInput );
	logInput->Update();

	// Subtract current bias field
	typename SubtractFilter::Pointer subFilter = SubtractFilter::New();
	subFilter->SetInput1( logInput->GetOutput() );
	subFilter->SetInput2( logModel->GetOutput() );
	subFilter->Update();

	// Estimate log bias function
	typename FieldEstimatorFilter::Pointer biasEstimator = FieldEstimatorFilter::New();
	typename FieldEstimatorFilter::ArrayType n;
	n.Fill( 5 );
	biasEstimator->SetNumberOfControlPoints( n );
	biasEstimator->SetInput( subFilter->GetOutput() );
	if ( m_MaskImage.IsNotNull() )
		biasEstimator->SetMaskImage( m_MaskImage );
	else {
		biasEstimator->SetMaskImage( modelFilter->GetEstimatedBackground() );
	}
	biasEstimator->Update();

	// Correct Input Image
	typename SubtractFilter::Pointer correct = SubtractFilter::New();
	correct->SetInput1( logInput->GetOutput() );
	correct->SetInput2( biasEstimator->GetOutput() );
	correct->Update();

	// Set up corrected image
	typename ExpFilter::Pointer exp = ExpFilter::New();
	exp->SetInput( correct->GetOutput() );
	exp->Update();
	m_CorrectedInput = exp->GetOutput();

	// Set new reference to sample
	m_Sample->SetImage( m_CorrectedInput );

	// Accumulate Bias
	typename AddFilter::Pointer biasAdder = AddFilter::New();
	biasAdder->SetInput1( m_BiasLog );
	biasAdder->SetInput2( biasEstimator->GetOutput() );
	biasAdder->Update();
	m_BiasLog = biasAdder->GetOutput();

	// Natural Units Bias
	typename ExpFilter::Pointer expBias = ExpFilter::New();
	expBias->SetInput( m_BiasLog );
	expBias->Update();
	m_CurrentBias = expBias->GetOutput();
}

template <class TInputVectorImage, class TProbabilityPixelType>
void ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::GenerateData() {
	int iteration = 0;
	m_CurrentIteration = 0;
	while (iteration < m_MaxIteration) {
		m_CurrentIteration = iteration;

		// E-step
		this->CalculateDensities();


		if( m_UseBiasCorrection ) // Bias field estimation
			this->UpdateBiasFieldEstimate();

		// M-step
		if (this->UpdateComponentParameters()) {
			this->UpdateProportions();
		} else {
			m_TerminationCode = CONVERGED;
			break;
		}


		// TODO Throw iteration event
		std::cout << "Iteration " << iteration << std::endl;
		for (unsigned int i = 0; i < m_ComponentVector.size(); i++) {
			std::cout << "\tClass[" << i << "," << m_Proportions[i] << "]: " << m_ComponentVector[i]->GetFullParameters() << std::endl;

		}

		++iteration;
	}

	m_TerminationCode = NOT_CONVERGED;
}

template <class TInputVectorImage, class TProbabilityPixelType>
const typename ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::MembershipFunctionVectorObjectType *
ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::GetOutput() const {
	size_t numberOfComponents = m_ComponentVector.size();
	MembershipFunctionVectorType & membershipFunctionsVector = m_MembershipFunctionsObject->Get();

	typename InputSampleType::MeasurementVectorSizeType measurementVectorSize = m_Sample->GetMeasurementVectorSize();

	typedef typename GaussianMembershipFunctionType::MeanVectorType MeanVectorType;

	MeanVectorType mean;
	itk::NumericTraits<MeanVectorType>::SetLength(mean, measurementVectorSize);

	typename GaussianMembershipFunctionType::CovarianceMatrixType covariance;
	covariance.SetSize(measurementVectorSize, measurementVectorSize);

	typename ComponentType::ParametersType parameters;

	for (size_t i = 0; i < numberOfComponents; ++i) {
		parameters = m_ComponentVector[i]->GetFullParameters();
		typename GaussianMembershipFunctionType::Pointer membershipFunction = GaussianMembershipFunctionType::New();
		membershipFunction->SetMeasurementVectorSize(measurementVectorSize);
		unsigned int parameterIndex = 0;
		for (unsigned int j = 0; j < measurementVectorSize; j++) {
			mean[j] = parameters[j];
			++parameterIndex;
		}

		for (unsigned int ii = 0; ii < measurementVectorSize; ++ii) {
			for (unsigned int jj = 0; jj < measurementVectorSize; ++jj) {
				covariance.GetVnlMatrix().put(ii, jj, parameters[parameterIndex]);
				++parameterIndex;
			}
		}

		membershipFunction->SetMean(mean);
		membershipFunction->SetCovariance(covariance);
		membershipFunctionsVector.push_back(membershipFunction.GetPointer());
	}

	return static_cast<const MembershipFunctionVectorObjectType *>(m_MembershipFunctionsObject);
}

template <class TInputVectorImage, class TProbabilityPixelType>
const typename ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::MembershipFunctionsWeightsArrayObjectType *
ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::GetMembershipFunctionsWeightsArray() const {
	size_t numberOfComponents = m_ComponentVector.size();
	ProportionVectorType & membershipFunctionsWeightVector = m_MembershipFunctionsWeightArrayObject->Get();

	membershipFunctionsWeightVector.SetSize(numberOfComponents);
	for (size_t i = 0; i < numberOfComponents; ++i) {
		membershipFunctionsWeightVector[i] = m_Proportions[i];
	}

	return static_cast<const MembershipFunctionsWeightsArrayObjectType *>(m_MembershipFunctionsWeightArrayObject);
}

template <class TInputVectorImage, class TProbabilityPixelType>
void ImageBiasedGMModelEMEstimator<TInputVectorImage,TProbabilityPixelType>::Update() {
	unsigned int numberOfComponents = m_ComponentVector.size();


	// Initialize a copy of input
	typedef itk::ImageDuplicator< InputVectorImageType > Duplicator;
	typename Duplicator::Pointer d = Duplicator::New();
	d->SetInputImage( m_Input );
	d->Update();
	m_CorrectedInput = d->GetOutput();

	// Initialize Sample
	m_Sample = InputSampleType::New();
	//m_Sample->SetImage( m_CorrectedInput );
	m_Sample->SetImage( m_Input );

	// Prepare mask for filters
	if( m_MaskImage.IsNotNull() ) {
		m_Sample->SetMask( m_MaskImage );
	}

	m_Sample->Initialize();

	// Initialize log bias field
	m_BiasLog = InputVectorImageType::New();
	m_BiasLog->SetRegions( m_Input->GetRequestedRegion() );
	m_BiasLog->CopyInformation( m_Input );
	m_BiasLog->SetNumberOfComponentsPerPixel( m_Input->GetNumberOfComponentsPerPixel() );
	m_BiasLog->Allocate();
	MeasurementVectorType m ( m_Input->GetNumberOfComponentsPerPixel() );
	m.Fill(0.0);
	m_BiasLog->FillBuffer( m );
	m_CurrentBias = InputVectorImageType::New();
	m_CurrentBias->SetRegions( m_Input->GetRequestedRegion() );
	m_CurrentBias->CopyInformation( m_Input );
	m_CurrentBias->SetNumberOfComponentsPerPixel( m_Input->GetNumberOfComponentsPerPixel() );
	m_CurrentBias->Allocate();
	m_CurrentBias->FillBuffer( m );


	// Initialize Model
	m_Proportions = m_InitialProportions;
	m_Posteriors.resize( numberOfComponents );

	for( size_t c = 0; c<numberOfComponents; c++) {
		m_Posteriors[c] = ProbabilityImageType::New();
		m_Posteriors[c]->SetRegions( m_Input->GetRequestedRegion() );
		m_Posteriors[c]->CopyInformation( m_Input );
		m_Posteriors[c]->Allocate();
		m_Posteriors[c]->FillBuffer(0.0);
		m_Posteriors[c]->Update();
	}

	this->GenerateData();
}
} // end of namespace Statistics
} // end of namespace itk


#endif /* IMAGEBIASEDGMMODELEMESTIMATOR_HXX_ */
