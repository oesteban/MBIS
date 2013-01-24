// --------------------------------------------------------------------------------------
// File:          MultispectralMLClassifierFilter.h
// Date:          Jan 24, 2013
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2013, code@oscaresteban.es (Oscar Esteban)
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

#ifndef MULTISPECTRALMLCLASSIFIERFILTER_H_
#define MULTISPECTRALMLCLASSIFIERFILTER_H_


#include <vector>

#include <itkImage.h>
#include <itkListSample.h>
#include <Numerics/Statistics/include/MaskedImageToListSampleAdaptor.h>
#include <itkDistanceToCentroidMembershipFunction.h>
#include <itkImageKmeansModelEstimator.h>

#include <itkNumericTraits.h>

#include "WeightedCentroidKdTreeImageGenerator.h"
#include <itkKdTree.h>
#include <itkKdTreeBasedKmeansEstimator.h>

#include <itkMinimumDecisionRule.h>
#include <itkSampleClassifierFilter.h>

#include <itkImageClassifierFilter.h>


namespace itk
{

template <class TInputComponent, class TMaskImage = TInputComponent >
class ITK_EXPORT MultispectralMLClassifierFilter:
	public MultispectralFilter< TInputComponent, itk::Image< unsigned char, TInputComponent::ImageDimension> >
{
public:
	/** Standard class typedefs */
	typedef MultispectralMLClassifierFilter                                              Self;
	typedef SmartPointer<Self>                                                     Pointer;
	typedef SmartPointer<const Self>                                               ConstPointer;

	typedef typename TInputComponent::PixelType                                    InputComponentPixelType;

	typedef Image< unsigned char, TInputComponent::ImageDimension>                 OutputImageType;
	typedef typename OutputImageType::Pointer                                      OutputImagePointer;
	typedef typename itk::ImageRegionIterator< OutputImageType >                   OutputImageIterator;

	typedef MultispectralFilter< TInputComponent, OutputImageType >                Superclass;
	typedef typename Superclass::InputVectorImageType                              InputVectorImageType;
	typedef typename InputVectorImageType::PixelType                               InputVectorPixelType;

    /** Mask Image typedefs */
    typedef TMaskImage                                                             MaskImageType;
    typedef typename MaskImageType::Pointer                                        MaskImagePointer;
    typedef typename MaskImageType::ConstPointer                                   MaskImageConstPointer;
    typedef typename MaskImageType::PixelType                                      MaskPixelType;
    typedef itk::Statistics::MaskedImageToListSampleAdaptor<MaskImageType, MaskImageType> MaskSampleType;
    typedef itk::Statistics::MaskedImageToListSampleAdaptor
    	                        < InputVectorImageType,MaskImageType >             InputSampleType;
    typedef typename InputSampleType::Pointer                                      InputSamplePointer;
    typedef typename InputSampleType::ConstPointer                                 InputSampleConstPointer;
    typedef typename InputSampleType::MeasurementVectorType                        MeasurementVectorType;

    /** Type used for representing the Mean values */
    typedef typename NumericTraits< InputComponentPixelType >::RealType            RealPixelType;

    typedef typename Superclass::ParametersType                                    ParametersType;
    typedef typename Superclass::ParametersVectorType                              ParametersVectorType;

    typedef itk::Statistics::DistanceToCentroidMembershipFunction
    		                           < MeasurementVectorType >                   MembershipFunctionType;
    typedef typename MembershipFunctionType::Pointer                               MembershipFunctionPointer;
    typedef itk::Statistics::MinimumDecisionRule                                   DecisionRuleType;
    typedef typename itk::Statistics::ImageClassifierFilter
    		< InputSampleType, InputVectorImageType, OutputImageType >             ImageClassifierType;
    typedef typename ImageClassifierType::ClassLabelVectorObjectType               ImageClassLabelVectorObjectType;
    typedef typename ImageClassifierType::ClassLabelVectorType                     ImageClassLabelVectorType;
    typedef typename ImageClassifierType::MembershipFunctionVectorObjectType       ImageMembershipFunctionVectorObjectType;
    typedef typename ImageClassifierType::MembershipFunctionVectorType             ImageMembershipFunctionVectorType;


	/** Method for creation through the object factory. */
	itkNewMacro(Self);
	/** Run-time type information (and related methods). */
	itkTypeMacro(MultispectralMLClassifierFilter, MultispectralFilter);

    itkSetObjectMacro( MaskImage, MaskImageType);
    itkGetObjectMacro( MaskImage, MaskImageType);

    itkGetConstMacro( Parameters, ParametersVectorType );
    void SetParameters( ParametersVectorType& p ) { m_Parameters = p; }

protected:
	MultispectralMLClassifierFilter();
	~MultispectralMLClassifierFilter() {};

	virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;

	/** Standard pipeline method. While this class does not implement a
	 * ThreadedGenerateData(), its GenerateData() delegates all
	 * calculations to an InputSampleFilter. Multithreading depends on
	 * that one  */
	void GenerateData();

	struct sortMeansClass {
	  bool operator () ( ParametersType i, ParametersType j) { return i[0]<j[0]; }
	} mySortMeans;

private:
	MultispectralMLClassifierFilter( const Self& ); //purposely not implemented
	void operator=(const Self& );       //purposely not implemented

	MaskImagePointer                          m_MaskImage;
	ParametersVectorType                m_Parameters;
};

}

#include "MultispectralMLClassifierFilter.hxx"



#endif /* MULTISPECTRALMLCLASSIFIERFILTER_H_ */
