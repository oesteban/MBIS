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

#ifndef __itkBiasedMixtureModelEMEstimator_h
#define __itkBiasedMixtureModelEMEstimator_h

#include "itkSimpleDataObjectDecorator.h"
#include "itkMixtureModelComponentBase.h"
#include "GaussianMembershipFunction.h"

using namespace itk;

namespace mfbs {
namespace Statistics {
/** \class BiasedMixtureModelEMEstimator
 *  \brief This class generates the parameter estimates for a mixture
 *  model using expectation maximization strategy.
 *
 * The first template argument is the type of the target sample
 * data. This estimator expects one or more mixture model component
 * objects of the classes derived from the
 * MixtureModelComponentBase. The actual component (or module)
 * parameters are updated by each component. Users can think this
 * class as a strategy or a integration point for the EM
 * procedure. The initial proportion (SetInitialProportions), the
 * input sample (SetSample), the mixture model components
 * (AddComponent), and the maximum iteration (SetMaximumIteration) are
 * required. The EM procedure terminates when the current iteration
 * reaches the maximum iteration or the model parameters converge.
 *
 * <b>Recent API changes:</b>
 * The static const macro to get the length of a measurement vector,
 * \c MeasurementVectorSize  has been removed to allow the length of a measurement
 * vector to be specified at run time. It is now obtained at run time from the
 * sample set as input. Please use the function
 * GetMeasurementVectorSize() to get the length.
 *
 * \sa MixtureModelComponentBase, GaussianMixtureModelComponent
 * \ingroup ITKStatistics
 *
 * \wiki
 * \wikiexample{Statistics/BiasedMixtureModelEMEstimator_2D,2D Gaussian Mixture Model Expectation Maximization}
 * \endwiki
 */

template<class TSample>
class BiasedMixtureModelEMEstimator: public itk::Object {
public:
	/** Standard class typedef */
	typedef BiasedMixtureModelEMEstimator Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	/** Standard macros */
	itkTypeMacro(BiasedMixtureModelEMEstimator, itk::Object);
	itkNewMacro (Self);

	/** TSample template argument related typedefs */
	typedef TSample SampleType;
	typedef typename TSample::MeasurementType MeasurementType;
	typedef typename TSample::MeasurementVectorType MeasurementVectorType;

	/** Typedef requried to generate dataobject decorated output that can
	 * be plugged into SampleClassifierFilter */
	typedef GaussianMembershipFunction<MeasurementVectorType> GaussianMembershipFunctionType;
	typedef typename GaussianMembershipFunctionType::Pointer GaussianMembershipFunctionPointer;

	typedef itk::Statistics::MembershipFunctionBase<MeasurementVectorType> MembershipFunctionType;
	typedef typename MembershipFunctionType::ConstPointer MembershipFunctionPointer;
	typedef std::vector<MembershipFunctionPointer> MembershipFunctionVectorType;
	typedef itk::SimpleDataObjectDecorator<MembershipFunctionVectorType> MembershipFunctionVectorObjectType;
	typedef typename MembershipFunctionVectorObjectType::Pointer MembershipFunctionVectorObjectPointer;

	/** Type of the mixture model component base class */
	typedef itk::Statistics::MixtureModelComponentBase<TSample> ComponentType;

	/** Type of the component pointer storage */
	typedef std::vector<ComponentType *> ComponentVectorType;

	/** Type of the membership function base class */
	typedef itk::Statistics::MembershipFunctionBase<MeasurementVectorType> ComponentMembershipFunctionType;

	/** Type of the array of the proportion values */
	typedef itk::Array<double> ProportionVectorType;

	/** Sets the target data that will be classified by this */
	void SetSample(const TSample *sample);

	/** Returns the target data */
	const TSample * GetSample() const;

	/** Set/Gets the initial proportion values. The size of proportion
	 * vector should be same as the number of component (or classes) */
	void SetInitialProportions(ProportionVectorType & propotion);

	const ProportionVectorType & GetInitialProportions() const;

	/** Gets the result proportion values */
	const ProportionVectorType & GetProportions() const;

	/** typedef for decorated array of proportion */
	typedef itk::SimpleDataObjectDecorator<ProportionVectorType> MembershipFunctionsWeightsArrayObjectType;
	typedef typename MembershipFunctionsWeightsArrayObjectType::Pointer MembershipFunctionsWeightsArrayPointer;

	/** Get method for data decorated Membership functions weights array */
	const MembershipFunctionsWeightsArrayObjectType * GetMembershipFunctionsWeightsArray() const;

	/** Set/Gets the maximum number of iterations. When the optimization
	 * process reaches the maximum number of interations, even if the
	 * class parameters aren't converged, the optimization process
	 * stops. */
	void SetMaximumIteration(int numberOfIterations);

	int GetMaximumIteration() const;

	/** Gets the current iteration. */
	int GetCurrentIteration() {
		return m_CurrentIteration;
	}

	/** Adds a new component (or class). */
	int AddComponent(ComponentType *component);

	/** Gets the total number of classes currently plugged in. */
	unsigned int GetNumberOfComponents() const;

	/** Runs the optimization process. */
	void Update();

	/** Termination status after running optimization */
	enum TERMINATION_CODE {
		CONVERGED = 0, NOT_CONVERGED = 1
	};

	/** Gets the termination status */
	TERMINATION_CODE GetTerminationCode() const;

	/** Gets the membership function specified by componentIndex
	 argument. */
	ComponentMembershipFunctionType * GetComponentMembershipFunction(int componentIndex) const;

	/** Output Membership function vector containing the membership functions with
	 * the final optimized parameters */
	const MembershipFunctionVectorObjectType * GetOutput() const;

protected:
	BiasedMixtureModelEMEstimator();
	virtual ~BiasedMixtureModelEMEstimator() {
	}
	void PrintSelf(std::ostream & os, Indent indent) const;

	bool CalculateDensities();

	double CalculateExpectation() const;

	bool UpdateComponentParameters();

	bool UpdateProportions();

	/** Starts the estimation process */
	void GenerateData();

private:
	/** Target data sample pointer*/
	const TSample *m_Sample;

	int m_MaxIteration;
	int m_CurrentIteration;

	TERMINATION_CODE m_TerminationCode;
	ComponentVectorType m_ComponentVector;
	ProportionVectorType m_InitialProportions;
	ProportionVectorType m_Proportions;
	std::vector<std::vector<double> > m_Ptly;

	MembershipFunctionVectorObjectPointer m_MembershipFunctionsObject;
	MembershipFunctionsWeightsArrayPointer m_MembershipFunctionsWeightArrayObject;
};
// end of class
}// end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "BiasedMixtureModelEMEstimator.hxx"
#endif

#endif
