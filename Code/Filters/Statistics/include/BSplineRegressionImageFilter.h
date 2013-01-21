// --------------------------------------------------------------------------------------
// File:          BSplineRegressionImageFilter.h
// Date:          May 18, 2012
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

#ifndef BSPLINEREGRESSIONIMAGEFILTER_H_
#define BSPLINEREGRESSIONIMAGEFILTER_H_

#include <itkVector.h>
#include <itkPointSet.h>
#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkBSplineScatteredDataPointSetToImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkComposeImageFilter.h>


using namespace itk;

namespace mfbs {
/** \class BSplineRegressionImageFilter
 *  \brief
 */

template< class TInputVectorImage >
class BSplineRegressionImageFilter: public ImageToImageFilter< TInputVectorImage, TInputVectorImage > {
public:
	typedef BSplineRegressionImageFilter                                     Self;
	typedef ImageToImageFilter< TInputVectorImage, TInputVectorImage >  Superclass;
	typedef SmartPointer< Self >                                             Pointer;
	typedef SmartPointer< const Self >                                      ConstPointer;

	itkTypeMacro(BSplineRegressionImageFilter, ImageToImageFilter);
	itkNewMacro(Self);

	itkStaticConstMacro( ImageDimension, unsigned int,  TInputVectorImage::ImageDimension );

	typedef TInputVectorImage                                                InputVectorImageType;
	typedef typename InputVectorImageType::Pointer                            InputVectorImagePointer;
	typedef typename InputVectorImageType::ConstPointer                       InputVectorImageConstPointer;
	typedef typename InputVectorImageType::PixelType                          MeasurementVectorType;
	typedef typename MeasurementVectorType::ValueType                         MeasurementType;

	typedef Image< MeasurementType, ImageDimension >                           ChannelImageType;
	typedef Image< float, ImageDimension >                                    MaskImageType;

	typedef VectorIndexSelectionCastImageFilter<InputVectorImageType, ChannelImageType>
	                                                                             InputSelectorType;
	typedef ComposeImageFilter< ChannelImageType >                             ChannelCombinatorType;

    /** B-spline smoothing filter argument typedefs */
	typedef Vector<MeasurementType, 1>                                          ScalarType;
	typedef PointSet<ScalarType, itkGetStaticConstMacro( ImageDimension )>      PointSetType;
	typedef Image<ScalarType, itkGetStaticConstMacro( ImageDimension )>         ScalarImageType;
	typedef typename PointSetType::Pointer                                     PointSetPointer;
	typedef typename PointSetType::PointType                                   PointType;

	/** B-sline filter typedefs */
	typedef BSplineScatteredDataPointSetToImageFilter<
	  PointSetType, ScalarImageType>                            BSplineFilterType;
	typedef typename BSplineFilterType::PointDataImageType    BiasFieldControlPointLatticeType;
	typedef typename BSplineFilterType::ArrayType             ArrayType;

	itkSetConstObjectMacro( MaskImage, MaskImageType );
	itkGetConstObjectMacro( MaskImage, MaskImageType );

	itkSetClampMacro( SplineOrder, size_t, 2, 10 );
	itkGetConstMacro( SplineOrder, size_t );

	itkSetMacro( NumberOfControlPoints, ArrayType );
	itkGetMacro( NumberOfControlPoints, ArrayType );

protected:
	BSplineRegressionImageFilter();
	virtual ~BSplineRegressionImageFilter(){}

	void PrintSelf(std::ostream & os, itk::Indent indent) const;

	/** Starts the estimation process */
	void GenerateData();
private:
	BSplineRegressionImageFilter(const Self &); //purposely not implemented
	void operator=(const Self &);                 //purposely not implemented

	InputVectorImageConstPointer m_Input;
	typename MaskImageType::ConstPointer m_MaskImage;

	// B-spline fitting parameters
	size_t       m_SplineOrder;
	ArrayType    m_NumberOfControlPoints;

}; // End of class
} // End namespace mfbs


#ifndef ITK_MANUAL_INSTANTIATION
#include "BSplineRegressionImageFilter.hxx"
#endif

#endif /* BSPLINEREGRESSIONIMAGEFILTER_H_ */
