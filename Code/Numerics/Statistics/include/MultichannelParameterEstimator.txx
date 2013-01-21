// --------------------------------------------------------------------------------------
// File:          MultichannelParameterEstimator.txx
// Date:          Jan 31, 2012
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

#ifndef MULTICHANNELPARAMETERESTIMATOR_TXX_
#define MULTICHANNELPARAMETERESTIMATOR_TXX_

namespace itk {

/**
 *
 */
template <class TInputComponent, class TProbabilityPixelType>
MultichannelParameterEstimator< TInputComponent, TProbabilityPixelType >
::MultichannelParameterEstimator() {
	// At least 1 inputs is necessary for a vector image.
	this->SetNumberOfRequiredInputs(0);
}

/**
 *
 */
template <class TInputComponent, class TProbabilityPixelType>
MultichannelParameterEstimator< TInputComponent, TProbabilityPixelType >
::~MultichannelParameterEstimator()
{
}


template <class TInputComponent, class TProbabilityPixelType>
void MultichannelParameterEstimator< TInputComponent, TProbabilityPixelType >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "ImageVector:" << std::endl;
  m_Input->Print( os, indent );
}

template <class TInputComponent, class TProbabilityPixelType>
void MultichannelParameterEstimator< TInputComponent, TProbabilityPixelType >
::GenerateData()
{
	Superclass::GenerateData();
	// Initialize sample
	m_Input = this->GetInputVectorImage();
    m_NumberOfComponents = m_Input->GetNumberOfComponentsPerPixel();
	m_Sample = InputSampleType::New();
	m_Sample->SetImage( m_Input );

	// Prepare mask for filters
	if( m_MaskImage.IsNotNull() ) {
		m_Sample->SetMask( m_MaskImage );
	}

	m_Sample->Initialize();

	// Compute Weighted Mean & Cov.
	// http://www.itk.org/Doxygen/html/classitk_1_1Statistics_1_1WeightedMeanSampleFilter.html
	// http://www.itk.org/Doxygen/html/classitk_1_1Statistics_1_1WeightedCovarianceSampleFilter.html
	// Weighted Mean & Covariance estimators
	typename CovarianceEstimatorType::Pointer covEst = CovarianceEstimatorType::New();
	covEst->SetInput(m_Sample);
	const typename CovarianceEstimatorType::MatrixDecoratedType * cov_out = covEst->GetCovarianceMatrixOutput();
	const typename CovarianceEstimatorType::MeasurementVectorDecoratedType * mean_out = covEst->GetMeanOutput();

	// Weights array generation for this class
	m_CachedWeights.SetSize( m_Sample->Size() );
	m_CachedWeights.Fill(1.0);
	typename InputSampleType::ConstIterator end = m_Sample->End();
	size_t pos = 0;
	for(typename InputSampleType::ConstIterator it = m_Sample->Begin(); it != end; ++it, pos++ ) {
		m_CachedWeights[pos] = *(m_Prior->GetBufferPointer() + it.GetInstanceIdentifier());
	}
	covEst->SetWeights( m_CachedWeights );

	try {
		covEst->Update();
	} catch( itk::ExceptionObject & ex ) {
		std::cerr << "Exception caught" << ex << std::endl;
	}


	// Generate Output
	m_OutputParameters = ParametersType( (m_NumberOfComponents + 1 ) * m_NumberOfComponents);
	for ( unsigned int k = 0; k < mean_out->Get().Size(); k++)
		m_OutputParameters[k] = mean_out->Get()[k];

	for ( unsigned int j = 0; j < cov_out->Get().Rows(); j++)
		for (unsigned int k = 0; k< cov_out->Get().Cols(); k++)
			m_OutputParameters[(j+1)*m_NumberOfComponents + k] = cov_out->Get()(k,j);
}


}


#endif /* MULTICHANNELPARAMETERESTIMATOR_TXX_ */
