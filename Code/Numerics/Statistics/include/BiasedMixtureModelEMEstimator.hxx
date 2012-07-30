// --------------------------------------------------------------------------------------
// File:          BiasedMixtureModelEMEstimator.hxx
// Date:          May 10, 2012
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

#ifndef __itkBiasedMixtureModelEMEstimator_hxx
#define __itkBiasedMixtureModelEMEstimator_hxx

#include "itkNumericTraits.h"
#include "BiasedMixtureModelEMEstimator.h"

namespace mfbs {
namespace Statistics {
template<class TSample>
BiasedMixtureModelEMEstimator<TSample>::BiasedMixtureModelEMEstimator() {
	m_TerminationCode = NOT_CONVERGED;

	m_MembershipFunctionsObject = MembershipFunctionVectorObjectType::New();
	m_MembershipFunctionsWeightArrayObject = MembershipFunctionsWeightsArrayObjectType::New();
	m_Sample = 0;
	m_MaxIteration = 100;
}

template<class TSample>
void BiasedMixtureModelEMEstimator<TSample>::PrintSelf(std::ostream & os, Indent indent) const {
	Superclass::PrintSelf(os, indent);
	os << indent << "Maximum Iteration: " << this->GetMaximumIteration() << std::endl;
	os << indent << "Sample: " << this->GetSample() << std::endl;
	os << indent << "Number Of Components: " << this->GetNumberOfComponents() << std::endl;
	for (unsigned int i = 0; i < this->GetNumberOfComponents(); i++) {
		os << indent << "Component Membership Function[" << i << "]: " << this->GetComponentMembershipFunction(i) << std::endl;
	}
	os << indent << "Termination Code: " << this->GetTerminationCode() << std::endl;
	os << indent << "Initial Proportions: " << this->GetInitialProportions() << std::endl;
	os << indent << "Proportions: " << this->GetProportions() << std::endl;
	os << indent << "Calculated Expectation: " << this->CalculateExpectation() << std::endl;
}

template<class TSample>
void BiasedMixtureModelEMEstimator<TSample>::SetMaximumIteration(int numberOfIterations) {
	m_MaxIteration = numberOfIterations;
}

template<class TSample>
int BiasedMixtureModelEMEstimator<TSample>::GetMaximumIteration() const {
	return m_MaxIteration;
}

template<class TSample>
void BiasedMixtureModelEMEstimator<TSample>::SetInitialProportions(ProportionVectorType & proportions) {
	m_InitialProportions = proportions;
}

template<class TSample>
const typename BiasedMixtureModelEMEstimator<TSample>::ProportionVectorType &
BiasedMixtureModelEMEstimator<TSample>::GetInitialProportions() const {
	return m_InitialProportions;
}

template<class TSample>
const typename BiasedMixtureModelEMEstimator<TSample>::ProportionVectorType &
BiasedMixtureModelEMEstimator<TSample>::GetProportions() const {
	return m_Proportions;
}

template<class TSample>
void BiasedMixtureModelEMEstimator<TSample>::SetSample(const TSample *sample) {
	m_Sample = sample;
}

template<class TSample>
const TSample *
BiasedMixtureModelEMEstimator<TSample>::GetSample() const {
	return m_Sample;
}

template<class TSample>
int BiasedMixtureModelEMEstimator<TSample>::AddComponent(ComponentType *component) {
	m_ComponentVector.push_back(component);
	return static_cast<int>(m_ComponentVector.size());
}

template<class TSample>
unsigned int BiasedMixtureModelEMEstimator<TSample>::GetNumberOfComponents() const {
	return m_ComponentVector.size();
}

template<class TSample>
typename BiasedMixtureModelEMEstimator<TSample>::TERMINATION_CODE BiasedMixtureModelEMEstimator<TSample>::GetTerminationCode() const {
	return m_TerminationCode;
}

template<class TSample>
typename BiasedMixtureModelEMEstimator<TSample>::ComponentMembershipFunctionType *
BiasedMixtureModelEMEstimator<TSample>::GetComponentMembershipFunction(int componentIndex) const {
	return (m_ComponentVector[componentIndex])->GetMembershipFunction();
}

template<class TSample>
bool BiasedMixtureModelEMEstimator<TSample>::CalculateDensities() {
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

	m_Ptly.resize(m_Sample->Size());

	typename TSample::ConstIterator last = m_Sample->End();

	typedef typename TSample::AbsoluteFrequencyType FrequencyType;
	FrequencyType frequency;
	FrequencyType zeroFrequency = NumericTraits < FrequencyType > ::Zero;
	typename TSample::MeasurementVectorType yVector;
	double density;
	double densitySum;
	double totalDensitySum;
	double minDouble = NumericTraits<double>::epsilon();

	SizeValueType measurementVectorIndex = 0;

	for (typename TSample::ConstIterator iter = m_Sample->Begin(); iter != last; ++iter, ++measurementVectorIndex) {
		m_Ptly[measurementVectorIndex].resize(numberOfComponents, 0.0);
		yVector = iter.GetMeasurementVector();
		densitySum = 0.0;

		// For every x belonging Lpm
		for (size_t x = 0; x < numberOfComponents; ++x) {
			double Px = m_Proportions[x];
			double Py_x = m_ComponentVector[x]->Evaluate(yVector);

			// This is P(x) * P(y=yVector|x, theta )
			density = Px * Py_x;
			tempWeights[x] = density;
			// This is the denominator of (25) in Bach2005
			densitySum += density;
		}

		for (size_t x = 0; x < numberOfComponents; ++x) {
			temp = tempWeights[x];

			// Final division to get P(k)(x|y_i,theta)
			if (densitySum > minDouble) // just to make sure temp does not blow up!
				temp /= densitySum;
			else
				temp = uniformWeight;

			// Save weight
			m_Ptly[measurementVectorIndex][x] = temp;
		}
	}

	return true;
}

template<class TSample>
double BiasedMixtureModelEMEstimator<TSample>::CalculateExpectation() const {
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

template<class TSample>
bool BiasedMixtureModelEMEstimator<TSample>::UpdateComponentParameters() {
	bool updated = false;
	size_t nComponents = m_ComponentVector.size();
	ComponentType *component;

	typename TSample::ConstIterator iter = m_Sample->Begin();
	typename TSample::ConstIterator last = m_Sample->End();

	size_t i = 0;
	double totalWeight;
	double temp;
	double minDouble = NumericTraits<double>::epsilon();

	while (iter != last) {
		totalWeight = 0.;

		for (size_t x = 0; x < nComponents; ++x) {
			totalWeight += m_Ptly[i][x];
		}

		for (size_t x = 0; x < nComponents; ++x) {
			temp = 0.;
			if (totalWeight > minDouble) // just to make sure temp does not blow up!
				temp = m_Ptly[i][x] / totalWeight;

			m_ComponentVector[x]->SetWeight(i, temp);
		}
		++iter;
		i++;
	}

	for (size_t x = 0; x < nComponents; ++x) {
		component = m_ComponentVector[x];
		component->Update();
		updated = updated || component->AreParametersModified();
	}

	return updated;
}

template<class TSample>
bool BiasedMixtureModelEMEstimator<TSample>::UpdateProportions() {
	// TODO Compute and cache the histogram so this is computed like SUM( h(bin) * GetWeight(j) );
	size_t numberOfComponents = m_ComponentVector.size();
	size_t sampleSize = m_Sample->Size();
	double totalFrequency = static_cast<double>(m_Sample->GetTotalFrequency());
	double totalSum = 0;
	bool updated = false;

	if (totalFrequency < NumericTraits<double>::epsilon()) {
		return updated;
	}

	double tempSum[numberOfComponents];

	for (size_t i = 0; i < numberOfComponents; ++i) {
		tempSum[i] = 0;
		for (size_t j = 0; j < totalFrequency; ++j) {
			tempSum[i] += (m_Ptly[j][i] * m_Sample->GetFrequency(j));
		}
		totalSum += tempSum[i];
	}

	for (size_t i = 0; i < numberOfComponents; ++i) {
		tempSum[i] /= totalSum;

		if (tempSum[i] != m_Proportions[i]) {
			m_Proportions[i] = tempSum[i];
			updated = true;
		}
	}

	return updated;
}

template<class TSample>
void BiasedMixtureModelEMEstimator<TSample>::GenerateData() {
	m_Proportions = m_InitialProportions;
	unsigned int nComponents = m_ComponentVector.size();

	int iteration = 0;
	m_CurrentIteration = 0;
	while (iteration < m_MaxIteration) {
		m_CurrentIteration = iteration;

		// E-step
		this->CalculateDensities();

		// M-step
		if (this->UpdateComponentParameters()) {
			this->UpdateProportions();
		} else {
			m_TerminationCode = CONVERGED;
			break;
		}

		std::cout << "Iteration " << iteration << std::endl;
		for (unsigned int i = 0; i < nComponents; i++) {
			std::cout << "\tClass[" << i << "," << m_Proportions[i] << "]: " << m_ComponentVector[i]->GetFullParameters() << std::endl;

		}

		++iteration;
	}

	m_TerminationCode = NOT_CONVERGED;
}

template<class TSample>
const typename BiasedMixtureModelEMEstimator<TSample>::MembershipFunctionVectorObjectType *
BiasedMixtureModelEMEstimator<TSample>::GetOutput() const {
	size_t numberOfComponents = m_ComponentVector.size();
	MembershipFunctionVectorType & membershipFunctionsVector = m_MembershipFunctionsObject->Get();

	typename SampleType::MeasurementVectorSizeType measurementVectorSize = m_Sample->GetMeasurementVectorSize();

	typename GaussianMembershipFunctionType::MeanVectorType mean;
	NumericTraits<typename GaussianMembershipFunctionType::MeanVectorType>::SetLength(mean, measurementVectorSize);

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

template<class TSample>
const typename BiasedMixtureModelEMEstimator<TSample>::MembershipFunctionsWeightsArrayObjectType *
BiasedMixtureModelEMEstimator<TSample>::GetMembershipFunctionsWeightsArray() const {
	size_t numberOfComponents = m_ComponentVector.size();
	ProportionVectorType & membershipFunctionsWeightVector = m_MembershipFunctionsWeightArrayObject->Get();

	membershipFunctionsWeightVector.SetSize(numberOfComponents);
	for (size_t i = 0; i < numberOfComponents; ++i) {
		membershipFunctionsWeightVector[i] = m_Proportions[i];
	}

	return static_cast<const MembershipFunctionsWeightsArrayObjectType *>(m_MembershipFunctionsWeightArrayObject);
}

template<class TSample>
void BiasedMixtureModelEMEstimator<TSample>::Update() {
	this->GenerateData();
}
} // end of namespace Statistics
} // end of namespace itk

#endif
