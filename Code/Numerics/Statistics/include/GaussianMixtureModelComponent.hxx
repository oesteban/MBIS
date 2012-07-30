// --------------------------------------------------------------------------------------
// File:          GaussianMixtureModelComponent.hxx
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


#ifndef __itkGaussianMixtureModelComponent_hxx
#define __itkGaussianMixtureModelComponent_hxx

#include <iostream>

#include "GaussianMixtureModelComponent.h"

namespace mfbs
{
namespace Statistics
{
template< class TSample >
GaussianMixtureModelComponent< TSample >
::GaussianMixtureModelComponent()
{
  m_MeanEstimator = MeanEstimatorType::New();
  m_CovarianceEstimator = CovarianceEstimatorType::New();
  m_GaussianMembershipFunction = NativeMembershipFunctionType::New();
  this->SetMembershipFunction( (MembershipFunctionType *)
                               m_GaussianMembershipFunction.GetPointer() );
  m_Mean.Fill(0.0);
  m_Covariance.SetIdentity();
}

template< class TSample >
void
GaussianMixtureModelComponent< TSample >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Mean: " << m_Mean << std::endl;
  os << indent << "Covariance: " << m_Covariance << std::endl;
  os << indent << "Mean Estimator: " << m_MeanEstimator << std::endl;
  os << indent << "Covariance Estimator: " << m_CovarianceEstimator << std::endl;
  os << indent << "GaussianMembershipFunction: " << m_GaussianMembershipFunction << std::endl;
}

template< class TSample >
void
GaussianMixtureModelComponent< TSample >
::SetSample(const TSample *sample)
{
  Superclass::SetSample(sample);

  m_MeanEstimator->SetInput(sample);
  m_CovarianceEstimator->SetInput(sample);

  const MeasurementVectorSizeType measurementVectorLength =
    sample->GetMeasurementVectorSize();
  m_GaussianMembershipFunction->SetMeasurementVectorSize(measurementVectorLength);

  NumericTraits<MeasurementVectorType>::SetLength(m_Mean, measurementVectorLength);
  m_Covariance.SetSize(measurementVectorLength, measurementVectorLength);

  m_Mean.Fill(NumericTraits< double >::Zero);

  m_Covariance.Fill(NumericTraits< double >::Zero);

  typename NativeMembershipFunctionType::MeanVectorType mean;

  NumericTraits<typename NativeMembershipFunctionType::MeanVectorType>::SetLength(mean,
    measurementVectorLength);

  for ( unsigned int i = 0; i < measurementVectorLength; ++i )
    {
    mean[i] = m_Mean[i];
    }

  m_GaussianMembershipFunction->SetMean(mean);
}

template< class TSample >
void
GaussianMixtureModelComponent< TSample >
::SetParameters(const ParametersType & parameters)
{
  Superclass::SetParameters(parameters);

  unsigned int paramIndex = 0;
  unsigned int i, j;

  bool changed = false;

  MeasurementVectorSizeType measurementVectorSize =
    this->GetSample()->GetMeasurementVectorSize();

  for ( i = 0; i < measurementVectorSize; i++ )
    {
    if ( m_Mean[i] != parameters[paramIndex] )
      {
      m_Mean[i] = parameters[paramIndex];
      changed = true;
      }

    ++paramIndex;
    }

  typename NativeMembershipFunctionType::MeanVectorType mean;

  NumericTraits<typename NativeMembershipFunctionType::MeanVectorType>::SetLength(mean,
    measurementVectorSize);

  for ( i = 0; i < measurementVectorSize; ++i )
    {
    mean[i] = m_Mean[i];
    }

  m_GaussianMembershipFunction->SetMean(mean);

  for ( i = 0; i < measurementVectorSize; i++ )
    {
    for ( j = 0; j < measurementVectorSize; j++ )
      {
      if ( m_Covariance.GetVnlMatrix().get(i, j) !=
           parameters[paramIndex] )
        {
        m_Covariance.GetVnlMatrix().put(i, j, parameters[paramIndex]);
        changed = true;
        }
      ++paramIndex;
      }
    }
  m_GaussianMembershipFunction->SetCovariance(m_Covariance);

  this->AreParametersModified(changed);
}

template< class TSample >
double
GaussianMixtureModelComponent< TSample >
::CalculateParametersChange()
{
  unsigned int i, j;

  typename MeanVectorType::MeasurementVectorType meanEstimate =
    m_MeanEstimator->GetMean();

  CovarianceMatrixType covEstimateDecoratedObject = m_CovarianceEstimator->GetOutput();
  typename CovarianceMatrixType::MeasurementVectorType covEstimate =  covEstimateDecoratedObject->et();

  double                    temp;
  double                    changes = 0.0;
  MeasurementVectorSizeType measurementVectorSize =
    this->GetSample()->GetMeasurementVectorSize();

  for ( i = 0; i < measurementVectorSize; i++ )
    {
    temp = m_Mean[i] - meanEstimate[i];
    changes += temp * temp;
    }

  for ( i = 0; i < measurementVectorSize; i++ )
    {
    for ( j = 0; j < measurementVectorSize; j++ )
      {
      temp = m_Covariance.GetVnlMatrix().get(i, j)
             - covEstimate.GetVnlMatrix().get(i, j);
      changes += temp * temp;
      }
    }

  changes = vcl_sqrt(changes);
  return changes;
}

template< class TSample >
void
GaussianMixtureModelComponent< TSample >
::GenerateData()
{
  MeasurementVectorSizeType measurementVectorSize =
    this->GetSample()->GetMeasurementVectorSize();

  this->AreParametersModified(false);

  const WeightArrayType & weights = this->GetWeights();

  typename TSample::ConstIterator iter = this->GetSample()->Begin();
  typename TSample::ConstIterator end =  this->GetSample()->End();

  typename TSample::MeasurementVectorType measurements;

  m_MeanEstimator->SetWeights(weights);
  m_MeanEstimator->Update();

  MeasurementVectorSizeType   i, j;
  double         temp;
  double         changes;
  bool           changed = false;
  ParametersType parameters = this->GetFullParameters();
  MeasurementVectorSizeType            paramIndex  = 0;

  typename MeanEstimatorType::MeasurementVectorType meanEstimate = m_MeanEstimator->GetMean();
  for ( i = 0; i < measurementVectorSize; i++ )
    {
    changes = vnl_math_abs( m_Mean[i] - meanEstimate[i] );

    if ( changes > this->GetMinimalParametersChange() )
      {
      changed = true;
      break;
      }
    }

  if ( changed )
    {
    m_Mean = meanEstimate;
    for ( paramIndex = 0; paramIndex < measurementVectorSize; paramIndex++ )
      {
      parameters[paramIndex] = meanEstimate[paramIndex];
      }
    this->AreParametersModified(true);
    }
  else
    {
    paramIndex = measurementVectorSize;
    }

  m_CovarianceEstimator->SetWeights(weights);
  m_CovarianceEstimator->Update();
  typename CovarianceEstimatorType::MatrixType covEstimate =
    m_CovarianceEstimator->GetCovarianceMatrix();

  changed = false;
  for ( i = 0; i < measurementVectorSize; i++ )
    {
    if( !changed )
      {
      for ( j = 0; j < measurementVectorSize; j++ )
        {
        temp = m_Covariance.GetVnlMatrix().get(i, j)
               - covEstimate.GetVnlMatrix().get(i, j);
        changes = vnl_math_abs( temp );
        if ( changes > this->GetMinimalParametersChange() )
          {
          changed = true;
          break;
          }
        }
      }
    }

  if ( changed )
    {
    m_Covariance = covEstimate;
    for ( i = 0; i < measurementVectorSize; i++ )
      {
      for ( j = 0; j < measurementVectorSize; j++ )
        {
        parameters[paramIndex] = covEstimate.GetVnlMatrix().get(i, j);
        ++paramIndex;
        }
      }
    this->AreParametersModified(true);
    }

  //THIS IS NEEDED TO update m_Mean and m_Covariance.SHOULD BE REMOVED
  paramIndex = 0;
  for ( i = 0; i < measurementVectorSize; i++ )
    {
    m_Mean[i] = parameters[paramIndex];
    ++paramIndex;
    }

  typename NativeMembershipFunctionType::MeanVectorType mean;
  NumericTraits<typename NativeMembershipFunctionType::MeanVectorType>::SetLength(mean,
    measurementVectorSize);

  for ( i = 0; i < measurementVectorSize; ++i )
    {
    mean[i] = m_Mean[i];
    }
  m_GaussianMembershipFunction->SetMean(mean);

  for ( i = 0; i < measurementVectorSize; i++ )
    {
    for ( j = 0; j < measurementVectorSize; j++ )
      {
      m_Covariance.GetVnlMatrix().put(i, j, parameters[paramIndex]);
      ++paramIndex;
      }
    }
  m_GaussianMembershipFunction->SetCovariance(m_Covariance);

  Superclass::SetParameters(parameters);
}
} // end of namespace Statistics
} // end of namespace itk

#endif
