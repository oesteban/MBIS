/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkMaskedImageToListSampleAdaptor_hxx
#define __itkMaskedImageToListSampleAdaptor_hxx

#include "MaskedImageToListSampleAdaptor.h"

namespace itk
{
namespace Statistics
{
template< class TImage, class TMask >
MaskedImageToListSampleAdaptor< TImage, TMask >
::MaskedImageToListSampleAdaptor()
{
  m_Image = 0;
  m_Mask = 0;
  m_TotalFrequency = NumericTraits< AbsoluteFrequencyType >::Zero;
  m_MeasurementVectorSize = 0;
}

template< class TImage, class TMask >
void MaskedImageToListSampleAdaptor< TImage, TMask >
::Initialize()
{
  if ( m_Image.IsNull() )
    {
    itkExceptionMacro("Image has not been set yet");
    }

  size_t imageSize = m_Image->GetLargestPossibleRegion().GetNumberOfPixels();
  m_MeasurementVectorSize = m_Image->GetNumberOfComponentsPerPixel();

  m_IdHolder.clear();

  if ( m_Mask.IsNotNull() ) {
	for ( size_t pos = 0; pos < imageSize; pos++) {
		if ( *(m_Mask->GetBufferPointer()+pos)>0.0) {
			m_IdHolder.push_back( pos );
		}
	}
  } else {
	for ( size_t pos = 0; pos < imageSize; pos++) {
		m_IdHolder.push_back( pos );
	}
  }
  m_TotalFrequency = m_IdHolder.size();

  this->Modified();
}

template< class TImage, class TMask >
inline const typename MaskedImageToListSampleAdaptor< TImage, TMask >::MeasurementVectorType &
MaskedImageToListSampleAdaptor< TImage, TMask >
::GetMeasurementVectorByIndex(size_t index) const
{

	MeasurementVectorType mv(m_MeasurementVectorSize);

	for( size_t i = 0; i < m_MeasurementVectorSize; i++) {
		mv[i] = (*m_Image->GetPixelContainer())[ m_IdHolder[index]*m_MeasurementVectorSize + i];
	}

	MeasurementVectorTraits::Assign( m_MeasurementVectorInternal, mv );

	return m_MeasurementVectorInternal;
}


template< class TImage, class TMask >
const typename MaskedImageToListSampleAdaptor< TImage, TMask >::MeasurementVectorType &
MaskedImageToListSampleAdaptor< TImage, TMask >
::GetMeasurementVector(InstanceIdentifier id) const
{
	MeasurementVectorType mv(m_MeasurementVectorSize);

//	for( size_t i = 0; i < m_MeasurementVectorSize; i++) {
//		mv[i] = (*m_Image->GetPixelContainer())[ m_IdHolder[id]*m_MeasurementVectorSize + i];
//	}

	for( size_t i = 0; i < m_MeasurementVectorSize; i++) {
		mv[i] = (*m_Image->GetPixelContainer())[ id * m_MeasurementVectorSize + i];
	}

	MeasurementVectorTraits::Assign( m_MeasurementVectorInternal, mv );

  return m_MeasurementVectorInternal;
}


template< class TImage, class TMask >
const typename MaskedImageToListSampleAdaptor< TImage, TMask >::MeasurementVectorType &
MaskedImageToListSampleAdaptor< TImage, TMask >
::GetRandomMeasurementVector() const {
	size_t u = rand() * RAND_MAX + rand();
	InstanceIdentifier id = ((u % m_TotalFrequency ) + m_TotalFrequency) % m_TotalFrequency;

	return GetMeasurementVector( m_IdHolder[id] );
}

template< class TImage, class TMask >
typename MaskedImageToListSampleAdaptor< TImage, TMask >::InstanceIdentifier
MaskedImageToListSampleAdaptor< TImage, TMask >
::Size() const
{
	 return m_TotalFrequency;
}


template< class TImage, class TMask >
inline typename MaskedImageToListSampleAdaptor< TImage, TMask >::AbsoluteFrequencyType
MaskedImageToListSampleAdaptor< TImage, TMask >
::GetFrequency(InstanceIdentifier id) const
{
  return NumericTraits< AbsoluteFrequencyType >::One;
}

template< class TImage, class TMask >
inline typename MaskedImageToListSampleAdaptor< TImage, TMask >::AbsoluteFrequencyType
MaskedImageToListSampleAdaptor< TImage, TMask >
::GetFrequencyByIndex(unsigned int index) const
{
  return NumericTraits< AbsoluteFrequencyType >::One;
}

template< class TImage, class TMask >
inline typename MaskedImageToListSampleAdaptor< TImage, TMask >::InstanceIdentifier
MaskedImageToListSampleAdaptor< TImage, TMask >
::GetInstanceIdentifier(unsigned int index) const
{
  if ( index >= m_IdHolder.size() )
    {
    itkExceptionMacro("Index out of range");
    }
  return m_IdHolder[index];
}

template< class TImage, class TMask >
inline void
MaskedImageToListSampleAdaptor< TImage, TMask >
::Swap(unsigned int index1, unsigned int index2)
{
  if ( index1 >= m_IdHolder.size()
       || index2 >= m_IdHolder.size() )
    {
    itkExceptionMacro("Index out of range");
    }

  InstanceIdentifier temp = m_IdHolder[index1];
  m_IdHolder[index1] = m_IdHolder[index2];
  m_IdHolder[index2] = temp;
  this->Modified();
}

template< class TImage, class TMask >
void
MaskedImageToListSampleAdaptor< TImage, TMask >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Image: ";
  if ( m_Image.IsNotNull() )
    {
    os << m_Image << std::endl;
    }
  else
    {
    os << "not set." << std::endl;
    }
  os << indent << "MeasurementVectorSize: ";
  if ( m_Image.IsNotNull() )
    {
    os << this->GetMeasurementVectorSize() << std::endl;
    }
  else
    {
    os << "not set." << std::endl;
    }
}

template< class TImage, class TMask >
void
MaskedImageToListSampleAdaptor< TImage, TMask >
::SetImage(const TImage *image)
{
  m_Image = image;
  this->Modified();
}

template< class TImage, class TMask >
void
MaskedImageToListSampleAdaptor< TImage, TMask >
::SetMask(const TMask *image)
{
  m_Mask = image;
  this->Modified();
}

template< class TImage, class TMask >
const TImage *
MaskedImageToListSampleAdaptor< TImage, TMask >
::GetImage() const
{
  if ( m_Image.IsNull() )
    {
    itkExceptionMacro("Image has not been set yet");
    }

  return m_Image.GetPointer();
}

template< class TImage, class TMask >
const TMask *
MaskedImageToListSampleAdaptor< TImage, TMask >
::GetMask() const
{
  if ( m_Mask.IsNull() )
    {
    itkExceptionMacro("Mask has not been set yet");
    }

  return m_Mask.GetPointer();
}

template< class TImage, class TMask >
typename MaskedImageToListSampleAdaptor< TImage, TMask >::TotalAbsoluteFrequencyType
MaskedImageToListSampleAdaptor< TImage, TMask >
::GetTotalFrequency() const
{
  return m_TotalFrequency;
}
} // end of namespace Statistics
} // end of namespace itk

#endif
