// --------------------------------------------------------------------------------------
// File:          MultichannelParameterEstimator.h
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

#ifndef MULTICHANNELPARAMETERESTIMATOR_H_
#define MULTICHANNELPARAMETERESTIMATOR_H_


#include <vector>
#include "Numerics/Statistics/include/MaskedImageToListSampleAdaptor.h"
#include <itkWeightedCovarianceSampleFilter.h>

#include "MultispectralFilter.h"

namespace itk
{

template <class TInputComponent, class TProbabilityPixelType = float >
class ITK_EXPORT MultichannelParameterEstimator:
	public MultispectralFilter< TInputComponent, itk::Image< unsigned char, TInputComponent::ImageDimension> >
{
public:
	/** Standard class typedefs */
	typedef MultichannelParameterEstimator                                                  Self;
	typedef SmartPointer<Self>                                                     Pointer;
	typedef SmartPointer<const Self>                                               ConstPointer;

	typedef typename TInputComponent::PixelType                                    InputComponentPixelType;

	typedef Image< unsigned char, TInputComponent::ImageDimension >                OutputImageType;
	typedef typename OutputImageType::Pointer                                      OutputImagePointer;
	typedef typename itk::ImageRegionIterator< OutputImageType >                   OutputImageIterator;

	typedef MultispectralFilter< TInputComponent, OutputImageType >                Superclass;

	typedef typename Superclass::InputVectorImageType                              InputVectorImageType;
	typedef typename itk::ImageRegionConstIteratorWithIndex
			                                      < InputVectorImageType >         InputImageIterator;

    typedef typename Superclass::ParametersType                                    ParametersType;

    /** Priors probability images typedefs */
    typedef TProbabilityPixelType                                                  ProbabilityPixelType;
    typedef itk::Image< ProbabilityPixelType, TInputComponent::ImageDimension >    ProbabilityImageType;
    typedef typename ProbabilityImageType::Pointer                                 ProbabilityImagePointer;
    typedef typename ProbabilityImageType::ConstPointer                            ProbabilityImageConstPointer;
    typedef typename itk::ImageRegionIterator< ProbabilityImageType >              ProbabilityImageIterator;
    typedef typename itk::ImageRegionConstIterator< ProbabilityImageType >         ProbabilityImageConstIterator;
/*
    typedef itk::Statistics::MaskedImageToListSampleAdaptor
    		            < ProbabilityImageType, ProbabilityImageType >             ProbabilitySampleType;
    typedef typename ProbabilitySampleType::Pointer                                ProbabilitySamplePointer;
    typedef typename ProbabilitySampleType::ConstPointer                           ProbabilitySampleConstPointer;*/

	typedef itk::Statistics::MaskedImageToListSampleAdaptor
	                    < InputVectorImageType, ProbabilityImageType >             InputSampleType;
	typedef typename InputSampleType::MeasurementVectorType                        MeasurementVectorType;

    /** Weighted statistical parameters estimation */
    /*typedef typename  itk::Statistics::WeightedMeanSampleFilter
    		                                    < InputSampleType >                MeanEstimatorType;*/
    typedef typename  itk::Statistics::WeightedCovarianceSampleFilter
    		                                    < InputSampleType >                CovarianceEstimatorType;
    typedef typename CovarianceEstimatorType::WeightArrayType                      WeightArrayType;
    typedef std::vector< WeightArrayType >                                         WeightArrayVectorType;


    /** Method to get the list sample, the generated output. Note that this does
     * not invoke Update(). You should have called update on this class to get
     * any meaningful output. */
    const InputSampleType * GetSample() const  { return m_Sample; }

	/** Method for creation through the object factory. */
	itkNewMacro(Self);
	/** Run-time type information (and related methods). */
	itkTypeMacro(MultichannelParameterEstimator, MultispectralFilter);


    itkSetConstObjectMacro( MaskImage, ProbabilityImageType);
    itkGetConstObjectMacro( MaskImage, ProbabilityImageType);

    itkSetConstObjectMacro( Prior, ProbabilityImageType);
    itkGetConstObjectMacro( Prior, ProbabilityImageType);

    itkGetMacro( CachedWeights, WeightArrayType);

    itkGetMacro( OutputParameters, ParametersType );

protected:
	MultichannelParameterEstimator();
	virtual ~MultichannelParameterEstimator();

	virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;

	  /** Standard pipeline method. While this class does not implement a
	   * ThreadedGenerateData(), its GenerateData() delegates all
	   * calculations to an InputSampleFilter. Multithreading depends on
	   * that one  */
	  void GenerateData();
private:
	MultichannelParameterEstimator( const Self& ); //purposely not implemented
	void operator=(const Self& );       //purposely not implemented


	typename InputVectorImageType::ConstPointer  m_Input;
	typename InputSampleType::Pointer            m_Sample;
	ProbabilityImageConstPointer                 m_MaskImage;
	ProbabilityImageConstPointer                 m_Prior;
	unsigned char                                m_NumberOfComponents;
	ParametersType                               m_OutputParameters;
	WeightArrayType                              m_CachedWeights;
};

}

#include "MultichannelParameterEstimator.txx"

#endif /* MULTICHANNELPARAMETERESTIMATOR_H_ */
