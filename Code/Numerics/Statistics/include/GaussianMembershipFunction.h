// --------------------------------------------------------------------------------------
// File:          itkBiasedGaussianMembershipFunction.h
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

#ifndef __itkGaussianMembershipFunction_h
#define __itkGaussianMembershipFunction_h

#include "itkMatrix.h"
#include "MembershipEnergyFunctionBase.h"

using namespace itk;

namespace mfbs {
namespace Statistics {
/** \class GaussianMembershipFunction
 * \brief GaussianMembershipFunction models class membership through a
 * multivariate Gaussian function.
 *
 * GaussianMembershipFunction is a subclass of MembershipFunctionBase
 * that models class membership (or likelihood) using a multivariate
 * Gaussian function. The mean and covariance structure of the
 * Gaussian are established using the methods SetMean() and
 * SetCovariance(). The mean is a vector-type that is the same
 * vector-type as the measurement vector but guarenteed to have a real
 * element type. For instance, if the measurement type is an
 * Vector<int,3>, then the mean is Vector<double,3>. If the
 * measurement type is a VariableLengthVector<float>, then the mean is
 * VariableLengthVector<double>. In contrast to this behavior, the
 * covariance is always a VariableSizeMatrix<double>.
 *
 * If the covariance is singular or nearly singular, the membership function
 * behaves somewhat like an impulse located at the mean. In this case,
 * we specify the covariance to be a diagonal matrix with large values
 * along the diagonal. This membership function, therefore,
 * will return small but differentiable values everywher and increase
 * sharply near the mean.
 *
 * \ingroup ITKStatistics
 */

template< class TMeasurementVector >
class GaussianMembershipFunction :	public MembershipEnergyFunctionBase< TMeasurementVector > {
public:
	/** Standard class typedefs */
	typedef GaussianMembershipFunction Self;
	typedef MembershipEnergyFunctionBase< TMeasurementVector > Superclass;
	typedef itk::SmartPointer< Self > Pointer;
	typedef itk::SmartPointer< const Self > ConstPointer;

	/** Standard macros */
	itkTypeMacro(GaussianMembershipFunction, MembershipEnergyFunctionBase);
	itkNewMacro(Self);

	/** SmartPointer class for superclass */
	typedef typename Superclass::Superclass::Pointer MembershipFunctionPointer;

	typedef typename Superclass::ParametersType ParametersType;

	/** Typedef alias for the measurement vectors */
	typedef TMeasurementVector MeasurementVectorType;

	/** Length of each measurement vector */
	typedef typename Superclass::MeasurementVectorSizeType MeasurementVectorSizeType;

	/** Type of the mean vector. RealType on a vector-type is the same
	 * vector-type but with a real element type.  */
	typedef typename itk::NumericTraits< MeasurementVectorType >::RealType MeasurementVectorRealType;
	typedef MeasurementVectorRealType MeanVectorType;

	/** Type of the covariance matrix */
	typedef itk::VariableSizeMatrix< double > CovarianceMatrixType;

	/** Set the mean of the Gaussian distribution. Mean is a vector type
	 * similar to the measurement type but with a real element type. */
	void SetMean(const MeanVectorType & mean);

	/** Get the mean of the Gaussian distribution. Mean is a vector type
	 * similar to the measurement type but with a real element type. */
	itkGetConstReferenceMacro(Mean, MeanVectorType);

	/** Set the covariance matrix. Covariance matrix is a
	 * VariableSizeMatrix of doubles. The inverse of the covariance
	 * matrix and the normlization term for the multivariate Gaussian
	 * are calculate whenever the covaraince matrix is changed. */
	void SetCovariance(const CovarianceMatrixType & cov);

	/* Get the covariance matrix. Covariance matrix is a
	 VariableSizeMatrix of doubles. */
	itkGetConstReferenceMacro(Covariance, CovarianceMatrixType);

	/* Get the inverse covariance matrix. Covariance matrix is a
	 VariableSizeMatrix of doubles. */
	itkGetConstReferenceMacro(InverseCovariance, CovarianceMatrixType);

	/** Set all parameters at once  	 */
	void SetParameters( const ParametersType& parameters );

	/** Evaluate the probability density of a measurement vector. */
	double Evaluate(const MeasurementVectorType & measurement) const;

	/** Evaluate the log-likelihood of a measurement vector. */
	double EvaluateEnergy(const MeasurementVectorType & measurement) const;

	/** Method to clone a membership function, i.e. create a new instance of
	 * the same type of membership function and configure its ivars to
	 * match. */
	virtual typename LightObject::Pointer InternalClone() const;

protected:
	GaussianMembershipFunction(void);
	virtual ~GaussianMembershipFunction(void) {}
	void PrintSelf(std::ostream & os, Indent indent) const;

	double ComputeMeasurement( const MeasurementVectorType & measurement ) const;

private:
	GaussianMembershipFunction(const Self &); //purposely not implemented
	void operator=(const Self &);//purposely not implemented

	MeanVectorType m_Mean;// mean
	CovarianceMatrixType m_Covariance;// covariance matrix

	// inverse covariance matrix. automatically calculated
	// when covariace matirx is set.
	CovarianceMatrixType m_InverseCovariance;

	// pre_factor (normalization term). automatically calculated
	// when covariace matirx is set.
	double m_PreFactor;
	double m_LogPreFactor;

	/** Boolean to cache whether the covarinace is singular or nearly singular */
	bool m_CovarianceNonsingular;
};
} // end of namespace Statistics
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "GaussianMembershipFunction.hxx"
#endif

#endif
