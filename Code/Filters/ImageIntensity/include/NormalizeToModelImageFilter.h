// --------------------------------------------------------------------------------------
// File:          NormalizeToModelImageFilter.h
// Date:          Jan 30, 2013
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

#ifndef NORMALIZETOMODELIMAGEFILTER_H_
#define NORMALIZETOMODELIMAGEFILTER_H_

#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkComposeImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>

using namespace itk;

namespace mfbs {

template<class TInputVectorImage >
class NormalizeToModelImageFilter: public ImageToImageFilter<TInputVectorImage, TInputVectorImage> {
public:
	/** Standard class typedefs */
	typedef NormalizeToModelImageFilter                              Self;
	typedef ImageToImageFilter<TInputVectorImage, TInputVectorImage> Superclass;
	typedef SmartPointer<Self>                                       Pointer;
	typedef SmartPointer<const Self>                                 ConstPointer;

	itkStaticConstMacro( Dimension, unsigned int, TInputVectorImage::ImageDimension );

	itkTypeMacro( NormalizeToModelImageFilter, ImageToImageFilter );
	itkNewMacro( NormalizeToModelImageFilter );

	typedef TInputVectorImage                                        InputVectorImageType;
	typedef typename InputVectorImageType::Pointer                   InputVectorImagePointer;
	typedef typename InputVectorImageType::ConstPointer              InputVectorImageConstPointer;
	typedef typename InputVectorImageType::PixelType                 MeasurementVectorType;
	typedef typename MeasurementVectorType::ValueType                MeasurementType;


	typedef InputVectorImageType                                     OutputImageType;
	typedef typename OutputImageType::Pointer                        OutputImagePointer;

	typedef itk::Image<MeasurementType, Dimension>                   InputChannelImageType;
	typedef itk::Image< float, Dimension >                      MaskImageType;

	typedef itk::VectorIndexSelectionCastImageFilter
			< InputVectorImageType, InputChannelImageType >          SelectionFilter;

	itkSetConstObjectMacro( MaskImage, MaskImageType );
	itkGetConstObjectMacro( MaskImage, MaskImageType );

	itkSetClampMacro(Lambda,float,0.0,1.0);
	itkGetConstMacro(Lambda, float );

	void SetReferenceImage(const InputVectorImageType* ref ) {
		this->SetNthInput( 0, const_cast< InputVectorImageType * >(ref) );
	}

	void SetNormalizeImage( const InputVectorImageType* norm ) {
		this->SetNthInput( 1, const_cast< InputVectorImageType * >(norm) );
	}

protected:
	NormalizeToModelImageFilter();
	~NormalizeToModelImageFilter() {};

	//virtual void PrintSelf(std::ostream& os, Indent indent) const;
	void GenerateData();
	//virtual void GenerateOutputInformation(void);

private:
	NormalizeToModelImageFilter(const Self&); //purposely not implemented
	void operator=(const Self&); //purposely not implemented

	typename MaskImageType::ConstPointer m_MaskImage;
	float m_Lambda;
};

}

#include "NormalizeToModelImageFilter.hxx"


#endif /* NORMALIZETOMODELIMAGEFILTER_H_ */
