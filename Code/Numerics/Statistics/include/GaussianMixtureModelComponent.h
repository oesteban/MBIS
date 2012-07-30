// --------------------------------------------------------------------------------------
// File:          GaussianMixtureModelComponent.h
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


#ifndef __itkGaussianMixtureModelComponent_h
#define __itkGaussianMixtureModelComponent_h

#include "itkMixtureModelComponentBase.h"
#include "GaussianMembershipFunction.h"
#include "WeightedCovarianceSampleFilter.h"

using namespace itk;

namespace mfbs
{
namespace Statistics
{
/** \class GaussianMixtureModelComponent
 * \brief is a component (derived from MixtureModelComponentBase) for
 * Gaussian class. This class is used in
 * ExpectationMaximizationMixtureModelEstimator.
 *
 * On every iteration of EM estimation, this class's GenerateData
 * method is called to compute the new distribution parameters.
 *
 * <b>Recent API changes:</b>
 * The static const macro to get the length of a measurement vector,
 * \c MeasurementVectorSize  has been removed to allow the length of a measurement
 * vector to be specified at run time. It is now obtained at run time from the
 * sample set as input. Please use the function
 * GetMeasurementVectorSize() to get the length.
 *
 * \sa MixtureModelComponentBase, ExpectationMaximizationMixtureModelEstimator
 * \ingroup ITKStatistics
 */

template< class TSample >
class GaussianMixtureModelComponent: public itk::Statistics::MixtureModelComponentBase< TSample > {
public:
  /**Standard class typedefs. */
  typedef GaussianMixtureModelComponent        Self;
  typedef itk::Statistics::MixtureModelComponentBase< TSample > Superclass;
  typedef itk::SmartPointer< Self >                 Pointer;
  typedef itk::SmartPointer< const Self >           ConstPointer;

  /**Standard Macros */
  itkTypeMacro(GaussianMixtureModelComponent, itk::Statistics::MixtureModelComponentBase);
  itkNewMacro(Self);

  /** Typedefs from the superclass */
  typedef typename Superclass::MeasurementVectorType     MeasurementVectorType;
  typedef typename Superclass::MeasurementVectorSizeType MeasurementVectorSizeType;
  typedef typename Superclass::MembershipFunctionType    MembershipFunctionType;
  typedef typename Superclass::WeightArrayType           WeightArrayType;
  typedef typename Superclass::ParametersType            ParametersType;

  /** Type of the membership function. Gaussian density function */
  typedef GaussianMembershipFunction< MeasurementVectorType >
  NativeMembershipFunctionType;

  /** Types of the mean and the covariance calculator that will update
   *  this component's distribution parameters */
  typedef itk::Statistics::WeightedMeanSampleFilter< TSample >       MeanEstimatorType;
  typedef itk::Statistics::WeightedCovarianceSampleFilter< TSample > CovarianceEstimatorType;

  /** Type of the mean vector */
  typedef typename MeanEstimatorType::OutputType MeanVectorType;

  /** Type of the covariance matrix */
  typedef typename CovarianceEstimatorType::OutputType CovarianceMatrixType;

  /** Sets the input sample */
  void SetSample(const TSample *sample);

  /** Sets the component's distribution parameters. */
  void SetParameters(const ParametersType & parameters);

protected:
  GaussianMixtureModelComponent();
  virtual ~GaussianMixtureModelComponent() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Returns the sum of squared changes in parameters between
   * iterations */
  double CalculateParametersChange();

  /** Computes the new distribution parameters */
  void GenerateData();

private:
  typename NativeMembershipFunctionType::Pointer m_GaussianMembershipFunction;

  typename MeanEstimatorType::MeasurementVectorType m_Mean;

  typename CovarianceEstimatorType::MatrixType m_Covariance;

  typename MeanEstimatorType::Pointer m_MeanEstimator;

  typename CovarianceEstimatorType::Pointer m_CovarianceEstimator;
};  // end of class
} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "GaussianMixtureModelComponent.hxx"
#endif

#endif
