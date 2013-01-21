// --------------------------------------------------------------------------------------
// File:    MultispectralFilter.h
// Date:    Sep 28, 2011
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

#ifndef MULTISPECTRALFILTER_H_
#define MULTISPECTRALFILTER_H_

#include <ostream>

#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkComposeImageFilter.h>

namespace itk {

/**
 * \class MultispectralFilter
 * \brief Interface class for multispectral segmentation and processing filters
 *
 * All input images are expected to have the same template parameters and have
 * the same size and origin.
 *
 * \sa Image
 * \sa VectorImage
 * \sa ImageToImageFilter
 * \sa ImageFeatureExtraction
 */

template <class TInputComponent, class TOutputImage >
class ITK_EXPORT MultispectralFilter: public itk::ImageToImageFilter< TInputComponent, TOutputImage >
{
public:
	/** Standard class typedefs */
	typedef MultispectralFilter                                  Self;
	typedef ImageToImageFilter< TInputComponent, TOutputImage >  Superclass;
	typedef SmartPointer<Self>                                   Pointer;
	typedef SmartPointer<const Self>                             ConstPointer;

	itkTypeMacro(MultispectralFilter, ImageToImageFilter);
	itkStaticConstMacro( Dimension, unsigned int,  TInputComponent::ImageDimension );
	itkStaticConstMacro( OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

	/** Image type information. */
	typedef TInputComponent                                        ComponentImageType;
	typedef TOutputImage                                           OutputImageType;

	/** Extract some information from the image types */
	typedef typename TOutputImage::PixelType                       OutputPixelType;
	typedef typename TOutputImage::Pointer                         OutputImagePointer;
	typedef typename TOutputImage::ConstPointer                    OutputImageConstPointer;

	typedef typename TInputComponent::PixelType                    ComponentPixelType;
	typedef typename TInputComponent::Pointer                      ComponentImagePointer;
	typedef typename TInputComponent::ConstPointer                 ComponentImageConstPointer;
	typedef typename std::vector< ComponentImageConstPointer >	   ComponentsVectorType;

	typedef ComposeImageFilter< ComponentImageType >               InputToVectorFilterType;
	typedef typename InputToVectorFilterType::Pointer              InputToVectorFilterPointer;
	typedef typename InputToVectorFilterType::OutputImageType      InputVectorImageType;
	typedef typename InputVectorImageType::PixelType               InputVectorPixelType;
	typedef typename InputVectorImageType::Pointer                 InputVectorImagePointer;
	typedef typename InputVectorImageType::ConstPointer            InputVectorImageConstPointer;

    typedef itk::Array< double >                                   ParametersType;
    typedef std::vector< ParametersType >                          ParametersVectorType;

	virtual void SetNthInput(unsigned int idx, const TInputComponent * image);
	virtual void SetInput( unsigned int idx, const TInputComponent * image);

	virtual void SetInputVector( const ComponentsVectorType& inputs );


	itkSetConstObjectMacro(InputVectorImage, InputVectorImageType);
	itkGetConstObjectMacro(InputVectorImage, InputVectorImageType);

protected:
	MultispectralFilter();
	virtual ~MultispectralFilter();

	//virtual void PrintSelf(std::ostream& os, Indent indent) const;

	//virtual void GenerateInputRequestedRegion();
	//virtual void GenerateOutputInformation(void);

	virtual void SetNthInput(unsigned int num, DataObject *input) {
		Superclass::SetNthInput(num, input);
	}

	void GenerateData();
private:
	MultispectralFilter( const Self& ); //purposely not implemented
	void operator=(const Self& );       //purposely not implemented

	InputToVectorFilterPointer          m_VectorFilter;
	InputVectorImageConstPointer        m_InputVectorImage;
};

}

#include "MultispectralFilter.txx"

#endif /* MULTISPECTRALFILTER_H_ */
