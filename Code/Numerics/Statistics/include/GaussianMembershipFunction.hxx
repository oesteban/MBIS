// --------------------------------------------------------------------------------------
// File:          itkOptGaussianMembershipFunction.hxx
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

#ifndef __itkGaussianMembershipFunction_hxx
#define __itkGaussianMembershipFunction_hxx

#include "GaussianMembershipFunction.h"

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_ldl_cholesky.h>

namespace mfbs {
namespace Statistics {
template<class TMeasurementVector>
GaussianMembershipFunction<TMeasurementVector>::GaussianMembershipFunction() {
	// Basic set-up for univariate by default
	NumericTraits<MeanVectorType>::SetLength(m_Mean, 1);
	m_Mean.Fill(0.0);
	m_PreFactor = 1.0 / vcl_sqrt(2.0 * vnl_math::pi); // default univariate
	m_LogPreFactor = 0.0;
	m_Covariance.SetSize(1, 1);
	m_Covariance.SetIdentity();
	m_InverseCovariance = m_Covariance;
	m_CovarianceNonsingular = true;
}

template<class TMeasurementVector>
void GaussianMembershipFunction<TMeasurementVector>::PrintSelf(std::ostream & os, Indent indent) const {
	Superclass::PrintSelf(os, indent);

	os << indent << "Mean: " << m_Mean << std::endl;
	os << indent << "Covariance: " << std::endl;
	os << m_Covariance.GetVnlMatrix();
	os << indent << "InverseCovariance: " << std::endl;
	os << indent << m_InverseCovariance.GetVnlMatrix();
	os << indent << "Prefactor: " << m_PreFactor << std::endl;
	os << indent << "Covariance nonsingular: " << (m_CovarianceNonsingular ? "true" : "false") << std::endl;
}

template<class TMeasurementVector>
void GaussianMembershipFunction<TMeasurementVector>::SetMean(const MeanVectorType & mean) {
	if (this->GetMeasurementVectorSize()) {
		itk::Statistics::MeasurementVectorTraits::Assert(mean, this->GetMeasurementVectorSize(),
				"GaussianMembershipFunction::SetMean(): Size of mean vector specified does not match the size of a measurement vector.");
	} else {
		// not already set, cache the size
		this->SetMeasurementVectorSize(mean.Size());
	}

	if (m_Mean != mean) {
		m_Mean = mean;
		this->Modified();
	}
}

template<class TMeasurementVector>
void GaussianMembershipFunction<TMeasurementVector>::SetCovariance(const CovarianceMatrixType & cov) {
	// Sanity check
	if (cov.GetVnlMatrix().rows() != cov.GetVnlMatrix().cols()) {
		itkExceptionMacro(<< "Covariance matrix must be square");
	}

	if ( this->GetMeasurementVectorSize() ){
		if ( cov.GetVnlMatrix().rows() != this->GetMeasurementVectorSize() ) {
		itkExceptionMacro(<< "Length of measurement vectors must be the same as the size of the covariance.");
		}
	} else {
		// not already set, cache the size
		this->SetMeasurementVectorSize( cov.GetVnlMatrix().rows() );
	}

	if (m_Covariance == cov) { // no need to copy the matrix, compute the inverse, or the normalization
		return;
	}

	// Compute diagonal and check that eigenvectors >= 0.0
	typedef vnl_diag_matrix<double>::iterator DiagonalIterator;
	typedef vnl_symmetric_eigensystem<double> Eigensystem;
	vnl_matrix< double > S = cov.GetVnlMatrix();
	Eigensystem* e = new Eigensystem( S );

	bool modified = false;
	DiagonalIterator itD = e->D.begin();
	while ( itD!= e->D.end() ) {
		if (*itD < 0) {
			*itD = 0.;
			modified = true;
		}
		itD++;
	}

	if (modified)
	m_Covariance = e->recompose();
	else
	m_Covariance = cov;

	delete e;

	// the inverse of the covariance matrix is first computed by SVD
	vnl_matrix_inverse< double > inv_cov( m_Covariance.GetVnlMatrix() );

	// the determinant is then costless this way
	double det = inv_cov.determinant_magnitude();

	if( det < 0.) {
		itkExceptionMacro( << "det (\\gamma) < 0" );
	}

	// FIXME Singurality Threshold for Covariance matrix: 1e-6 is an arbitrary value!!!
	const double singularThreshold = 1.0e-6;
	m_CovarianceNonsingular = ( det > singularThreshold );

	if( m_CovarianceNonsingular ) {
		// allocate the memory for m_InverseCovariance matrix
		m_InverseCovariance.GetVnlMatrix() = inv_cov.inverse();

	} else {
		// Perform cholesky diagonalization and select the semi-positive aproximation
		itkWarningMacro( "Warning: covariance is singular, setting diagonal covariance" );
		vnl_ldl_cholesky* chol = new vnl_ldl_cholesky( m_Covariance.GetVnlMatrix() );
		vnl_vector< double > D( chol->diagonal() );
		det = dot_product( D, D );
		vnl_matrix_inverse< double > R (chol->upper_triangle());
		m_InverseCovariance.GetVnlMatrix() = R.inverse();
		m_CovarianceNonsingular = true;
	}

	// calculate coefficient C of multivariate gaussian
	double d = static_cast< double >( this->GetMeasurementVectorSize() );
	m_PreFactor = 1.0 / vcl_pow( vcl_pow( 2.0 * vnl_math::pi, d ) * det, 0.5);
	m_LogPreFactor = vcl_log( det ) + d * vcl_log( 2.0*vnl_math::pi );

	this->Modified();
}

template<class TMeasurementVector>
void GaussianMembershipFunction<TMeasurementVector>::SetParameters(const ParametersType& parameters ) {
	MeanVectorType mean;
	CovarianceMatrixType cov;

	size_t measurementVectorSize = 0.5 * (vcl_sqrt( 1 + 4 * parameters.Size() ) - 1 );

	mean.SetSize( measurementVectorSize );
	cov.SetSize( measurementVectorSize, measurementVectorSize );

	for( size_t i = 0; i < measurementVectorSize; i++ ) {
		mean[i] = parameters[i];
	}

	SetMean( mean );

	for( size_t i = 0; i < measurementVectorSize; i++ ) {
		for ( size_t j = 0; j < measurementVectorSize; j++) {
			cov(i,j) = parameters[j*measurementVectorSize+ i + measurementVectorSize];
		}
	}

	SetCovariance( cov );
}

template<class TMeasurementVector>
inline double GaussianMembershipFunction<TMeasurementVector>
::ComputeMeasurement(const MeasurementVectorType & measurement) const {
	const MeasurementVectorSizeType measurementVectorSize = this->GetMeasurementVectorSize();
	double argument = 0.;

	// Compute ( y - mean )
	vnl_vector<double> tempVector(measurementVectorSize);

	for (MeasurementVectorSizeType i = 0; i < measurementVectorSize; ++i) {
		tempVector[i] = measurement[i] - m_Mean[i];
	}

	// temp = ( y - mean )^t * InverseCovariance * ( y - mean )
	argument = dot_product(tempVector, m_InverseCovariance.GetVnlMatrix() * tempVector);

	return argument;
}


template<class TMeasurementVector>
inline double GaussianMembershipFunction<TMeasurementVector>::Evaluate(const MeasurementVectorType & measurement) const {
	return m_PreFactor * vcl_exp( -0.5 * this->ComputeMeasurement( measurement ) );
}

template<class TMeasurementVector>
inline double GaussianMembershipFunction<TMeasurementVector>::EvaluateEnergy(const MeasurementVectorType & measurement) const {
	return 0.5 * ( m_LogPreFactor + this->ComputeMeasurement( measurement ) );
}

template<class TVector>
typename LightObject::Pointer GaussianMembershipFunction<TVector>
::InternalClone() const {
	LightObject::Pointer loPtr = Superclass::InternalClone();
	typename Self::Pointer membershipFunction = dynamic_cast<Self *>(loPtr.GetPointer());
	if (membershipFunction.IsNull()) {
		itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass() << " failed.");
	}

	membershipFunction->SetMeasurementVectorSize(this->GetMeasurementVectorSize());
	membershipFunction->SetMean(this->GetMean());
	membershipFunction->SetCovariance(this->GetCovariance());

	return loPtr;
}

} // end namespace Statistics
} // end of namespace itk

#endif
