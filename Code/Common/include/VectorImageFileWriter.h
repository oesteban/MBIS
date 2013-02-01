// --------------------------------------------------------------------------------------
// File:          VectorImageFileWriter.h
// Date:          Jan 31, 2013
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

#ifndef VECTORIMAGEFILEWRITER_H_
#define VECTORIMAGEFILEWRITER_H_

// Include headers
#include <itkObject.h>
#include <itkImageFileWriter.h>

// Namespace declaration

namespace mfbs {
/** \class VectorImageFileWriter
 *  \brief This class
 *
 *  Long description
 *
 *  \ingroup
 */

template< typename TVectorImage >
class VectorImageFileWriter: public itk::Object {
public:
	typedef VectorImageFileWriter   Self;
	typedef itk::Object                   Superclass;
	typedef itk::SmartPointer<Self>       Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro( VectorImageFileWriter, itk::Object );
	itkNewMacro( Self );

	itkStaticConstMacro( Dimension, unsigned int, TVectorImage::ImageDimension );

	typedef TVectorImage                                      VectorImageType;
	typedef typename VectorImageType::ConstPointer            VectorImagePointer;
	typedef typename VectorImageType::PixelType               VectorType;
	typedef typename VectorType::ValueType                    ValueType;
	typedef typename itk::Image< ValueType, Dimension >       ComponentType;

	typedef itk::VectorIndexSelectionCastImageFilter
			         < VectorImageType, ComponentType >       SelectionFilter;


	itkSetConstObjectMacro(Input,VectorImageType);
	itkSetStringMacro(FileName);

	void Update() const {
		typedef itk::Image<ValueType,Dimension+1> OutputImageType;
		typename OutputImageType::Pointer out = OutputImageType::New();
		typename OutputImageType::SizeType outSize;
		typename OutputImageType::SpacingType outSpacing;

		size_t nComponents = m_Input->GetNumberOfComponentsPerPixel();

		for( size_t i = 0; i<Dimension; i++) {
			outSize[i] = m_Input->GetLargestPossibleRegion().GetSize()[i];
			outSpacing[i] = m_Input->GetSpacing()[i];
		}

		outSize[Dimension] = nComponents;
		outSpacing[Dimension] = 1.0;
		out->SetRegions( outSize );
		out->SetSpacing( outSpacing );
		out->Allocate();
		out->FillBuffer(0.0);
		ValueType* buffer = out->GetBufferPointer();


		size_t nVect = m_Input->GetLargestPossibleRegion().GetNumberOfPixels();

		typename SelectionFilter::Pointer sel = SelectionFilter::New();
		sel->SetInput( m_Input );

		size_t pix = 0;
		for( size_t comp = 0; comp < nComponents; comp++ ) {
			sel->SetIndex( comp );
			sel->Update();
			typename ComponentType::Pointer im = sel->GetOutput();
			size_t nPix = im->GetLargestPossibleRegion().GetNumberOfPixels();

			const ValueType* compBuffer = im->GetBufferPointer();

			for( size_t i = 0; i<nPix; i++ ) {
				ValueType v = *(compBuffer+i);
				*(buffer+pix) = v;
				 pix++;
			}
		}

		typename itk::ImageFileWriter<OutputImageType>::Pointer w = itk::ImageFileWriter<OutputImageType>::New();
		w->SetInput( out );
		w->SetFileName( m_FileName );
		w->Update();
	}

protected:
	VectorImageFileWriter(){}
	~VectorImageFileWriter(){}

	void PrintSelf( std::ostream& os, itk::Indent indent) const {
		Superclass::Print(os, indent);
	}

private:
	VectorImageFileWriter( const Self &); // purposely not implemented
	void operator=(const Self &); // purposely not implemented

	std::string m_FileName;
	VectorImagePointer m_Input;
}; // End of class VectorImageFileWriter
} // End of namespace





#endif /* VECTORIMAGEFILEWRITER_H_ */
