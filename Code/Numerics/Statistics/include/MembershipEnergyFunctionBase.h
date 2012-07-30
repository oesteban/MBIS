/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __MembershipFunctionBase_h
#define __MembershipFunctionBase_h

#include "itkMembershipFunctionBase.h"
#include "itkArray.h"

using namespace itk::Statistics;

namespace mfbs {
namespace Statistics {
/** \class MembershipFunctionBase
 * \brief MembershipFunctionBase defines common interfaces
 * for membership functions.
 *
 * MembershipFunctionBase is a subclass of FunctionBase which
 * restricts the function type to be a membership function. Membership
 * functions provide a mapping from an arbitrary domain to a set of
 * real numbers. Membership functions are typically used to model or
 * approximate likelihood functions, \f$p( x | i )\f$, i.e. the
 * probability of the measurement \f$x\f$ belonging to a class
 * \f$i\f$.
 *
 * The Statistics framework models random variables \f$x\f$ as
 * vectors. Typical uses of MembershipFunctions include templating
 * over a FixedArray, Array, Vector, or VariableLengthVector.
 *
 * The Evaluate() method returns the membership rank or likelihood
 * that the measurement belongs to the class represented by this
 * membership function.
 *
 * Evaluations of a single measurement across of set MembershipFunctions
 * can then be passed to a DecisionRule in order to establish
 * class (or group) assignment.
 *
 * \ingroup ITKStatistics
 */

template< class TVector >
class ITK_EXPORT MembershipEnergyFunctionBase: public MembershipFunctionBase< TVector >
{
public:
	/** Standard class typedefs */
	typedef MembershipEnergyFunctionBase Self;
	typedef MembershipFunctionBase< TVector > Superclass;
	typedef itk::SmartPointer< Self > Pointer;
	typedef itk::SmartPointer< const Self > ConstPointer;

	/** Standard macros */
	itkTypeMacro(MembershipFunctionBase, MembershipFunctionBase);

	/** MeasurementVector typedef support */
	typedef TVector MeasurementVectorType;

	/** Typedef for the length of each measurement vector */
	typedef unsigned int MeasurementVectorSizeType;


	typedef itk::Array< double > ParametersType;

	/** Method to get energy score of an entity
	 * or measurement. EvaluateEnergy() maps from a vector measurement type
	 * to a real number. */
	virtual double EvaluateEnergy(const MeasurementVectorType & x) const = 0;

	virtual void SetParameters(const ParametersType& p ) = 0;

protected:
	MembershipEnergyFunctionBase() {
		m_MeasurementVectorSize = itk::NumericTraits<MeasurementVectorType>::GetLength(
				MeasurementVectorType() );
	}
	virtual ~MembershipEnergyFunctionBase(void) {}

private:
	MembershipEnergyFunctionBase(const Self &); //purposely not implemented
	void operator=(const Self &);//purposely not implemented

	 MeasurementVectorSizeType m_MeasurementVectorSize;
}; // end of class
} // end of namespace Statistics
} // end namespace mfbs

#endif
