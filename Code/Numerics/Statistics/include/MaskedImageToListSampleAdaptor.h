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
#ifndef __itkMaskedImageToListSampleAdaptor_h
#define __itkMaskedImageToListSampleAdaptor_h

#include <typeinfo>

#include "itkImage.h"
#include "itkPixelTraits.h"
#include "itkListSample.h"
#include "itkSubsample.h"
#include "itkSmartPointer.h"
#include "itkImageRegionIterator.h"
#include "itkMeasurementVectorTraits.h"

namespace itk {
namespace Statistics {
/** \class MaskedImageToListSampleAdaptor
 *  \brief This class provides ListSample interface to ITK Image
 *
 * After calling SetImage( const Image * ) method to plug in the image object,
 * users can use Sample interfaces to access Image data. The resulting data
 * are a list of measurement vectors.
 *
 * The measurment vector type is determined from the image pixel type. This class
 * handles images with scalar, fixed array or variable length vector pixel types.
 *
 * \sa Sample, ListSample
 * \ingroup ITKStatistics
 *
 * \wiki
 * \wikiexample{Statistics/MaskedImageToListSampleAdaptor,Create a list of samples from an image without duplicating the data}
 * \endwiki
 */

template< class TImage, class TMask >
class ITK_EXPORT MaskedImageToListSampleAdaptor:
public Subsample< ListSample<typename MeasurementVectorPixelTraits< typename TImage::PixelType >::MeasurementVectorType> >
{
public:
	/** Standard class typedefs */
	typedef MaskedImageToListSampleAdaptor Self;

	typedef ListSample<typename MeasurementVectorPixelTraits< typename TImage::PixelType >::MeasurementVectorType> SampleType;

	typedef Subsample< SampleType >	Superclass;

	typedef SmartPointer< Self > Pointer;
	typedef SmartPointer< const Self > ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro(MaskedImageToListSampleAdaptor, Subsample);

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Image typedefs */
	typedef TImage ImageType;
	typedef typename ImageType::Pointer ImagePointer;
	typedef typename ImageType::ConstPointer ImageConstPointer;
	typedef typename ImageType::IndexType IndexType;
	typedef typename ImageType::PixelType PixelType;
	typedef typename ImageType::PixelContainerConstPointer PixelContainerConstPointer;

	/** Image Iterator typedef support */

	/** Mask typedefs */
	typedef TMask MaskType;
	typedef typename MaskType::Pointer MaskPointer;
	typedef typename MaskType::ConstPointer MaskConstPointer;
	typedef typename MaskType::IndexType MaskIndexType;
	typedef typename MaskType::PixelType MaskPixelType;
	typedef typename MaskType::PixelContainerConstPointer MaskPixelContainerConstPointer;

	/** Mask Iterator typedef support */
	typedef ImageRegionIterator< MaskType > MaskIteratorType;
	typedef ImageRegionConstIterator< MaskType > MaskConstIteratorType;
	typedef PixelTraits< typename TMask::PixelType > MaskPixelTraitsType;

	/** Superclass typedefs for Measurement vector, measurement,
	 * Instance Identifier, frequency, size, size element value */
	typedef MeasurementVectorPixelTraits< PixelType > MeasurementPixelTraitsType;
	typedef typename MeasurementPixelTraitsType::MeasurementVectorType MeasurementVectorType;

	typedef MeasurementVectorTraitsTypes< MeasurementVectorType > MeasurementVectorTraitsType;
	typedef typename MeasurementVectorTraitsType::ValueType MeasurementType;

	typedef typename Superclass::AbsoluteFrequencyType AbsoluteFrequencyType;
	typedef typename Superclass::TotalAbsoluteFrequencyType TotalAbsoluteFrequencyType;
	typedef typename Superclass::MeasurementVectorSizeType MeasurementVectorSizeType;
	typedef typename Superclass::InstanceIdentifier InstanceIdentifier;

	typedef std::vector< size_t >                    InstanceIdentifierHolder;
	typedef InstanceIdentifierHolder::const_iterator OffsetsHolderIterator;
	typedef typename TImage::PixelType*              PixelBufferPointer;
	typedef typename TImage::PixelContainer*         PixelContainerPointer;

	typedef MeasurementVectorType ValueType;

	/** Method to set the image */
	void SetImage(const TImage *image);

	/** Method to get the image */
	const TImage * GetImage() const;

	/** Method to set the Mask */
	void SetMask(const TMask *image);

	/** Method to get the Mask */
	const TMask * GetMask() const;

	virtual const InstanceIdentifierHolder & GetIdHolder() const
	{ \
	    return this->m_IdHolder;
	}

	/** returns the number of measurement vectors in this container */
	InstanceIdentifier Size() const;

	InstanceIdentifier GetInstanceIdentifier( unsigned int index ) const;

	/** Initialize method */
	void Initialize();

	/** Initialize the subsample with all instances of the sample */
	void InitializeWithAllInstances() { this->Initialize(); }

	const MeasurementVectorType & GetMeasurementVectorByIndex(size_t index) const;

	inline const size_t GetOffsetByIndex(size_t index) const {
		return m_IdHolder[index];
	}

	const MeasurementVectorType & GetRandomMeasurementVector() const;

	/** method to return measurement vector for a specified id */
	const MeasurementVectorType & GetMeasurementVector(InstanceIdentifier id) const;


	MeasurementVectorSizeType GetMeasurementVectorSize() const
	{
		// some filter are expected that this method returns something even if the
		// input is not set. This won't be the right value for a variable length vector
		// but it's better than an exception.
		if( m_Image.IsNull() )
		{
			return Superclass::GetMeasurementVectorSize();
		}
		else
		{
			return m_MeasurementVectorSize;
		}
	}

	/** method to return frequency for a specified id */
	AbsoluteFrequencyType GetFrequency(InstanceIdentifier id) const;
	AbsoluteFrequencyType GetFrequencyByIndex(unsigned int index) const;

	void Swap(unsigned int index1, unsigned int index2);

	/** method to return the total frequency */
	TotalAbsoluteFrequencyType GetTotalFrequency() const;

	/** \class ConstIterator
	 *  \brief Const Iterator
	 * \ingroup ITKStatistics
	 */
	class ConstIterator
	{
		friend class MaskedImageToListSampleAdaptor;
	public:

		ConstIterator(const MaskedImageToListSampleAdaptor *adaptor)
		{
			*this = adaptor->Begin();
		}

		ConstIterator(const ConstIterator & iter)
		{
			m_Iter = iter.m_Iter;
			m_ImageBuffer = iter.m_ImageBuffer;
			m_PixelSize = iter.m_PixelSize;
		}

		ConstIterator & operator=(const ConstIterator & iter)
		{
			m_Iter = iter.m_Iter;
			m_ImageBuffer = iter.m_ImageBuffer;
			m_PixelSize = iter.m_PixelSize;
			return *this;
		}

		AbsoluteFrequencyType GetFrequency() const
		{
			return NumericTraits< AbsoluteFrequencyType >::One;
		}

		const MeasurementVectorType & GetMeasurementVector() const
		{
			MeasurementVectorType mv(this->m_PixelSize);
			for( size_t i = 0; i < this->m_PixelSize; i++) {
				mv[i] = (*m_ImageBuffer)[(*m_Iter)*this->m_PixelSize + i];
			}

			MeasurementVectorTraits::Assign( this->m_MeasurementVectorCache, mv );
			return this->m_MeasurementVectorCache;
		}

		inline size_t GetInstanceIdentifier() const
		{
			return *m_Iter;
		}

		inline ConstIterator & operator++()
		{
			m_Iter++;
			return *this;
		}

		bool operator<(const ConstIterator & it)
		{
			return ( *m_Iter < *(it.m_Iter) );
		}

		bool operator!=(const ConstIterator & it)
		{
			return ( m_Iter != it.m_Iter );
		}

		bool operator==(const ConstIterator & it)
		{
			return ( m_Iter == it.m_Iter );
		}

	protected:
		// This method should only be available to the ListSample class
		ConstIterator( OffsetsHolderIterator iter, PixelContainerPointer ptr, size_t size)
		{
			this->m_Iter = iter;
			this->m_ImageBuffer = ptr;
			this->m_PixelSize = size;
		}

		// This method is purposely not implemented
		ConstIterator();
	private:
		mutable MeasurementVectorType m_MeasurementVectorCache;
		PixelContainerPointer m_ImageBuffer;
		OffsetsHolderIterator m_Iter;
		size_t m_PixelSize;
	};

	/** \class Iterator
	 *  \brief Iterator
	 * \ingroup ITKStatistics
	 */
	class Iterator:public ConstIterator
	{
		friend class MaskedImageToListSampleAdaptor;
	public:

		Iterator(Self *adaptor):ConstIterator(adaptor)
		{}

		Iterator(const Iterator & iter):ConstIterator(iter)
		{}

		Iterator & operator=(const Iterator & iter)
		{
			this->ConstIterator::operator=(iter);
			return *this;
		}

	protected:
		// To ensure const-correctness these method must not be in the public API.
		// The are purposly not implemented, since they should never be called.
		Iterator();
		Iterator(const Self *adaptor);
		//Iterator( OffsetsHolderIterator iter );
		//Iterator(ImageConstIteratorType iter, InstanceIdentifier iid);
		Iterator(const ConstIterator & it);
		ConstIterator & operator=(const ConstIterator & it);

		Iterator(OffsetsHolderIterator iter, PixelContainerPointer ptr, size_t size ):ConstIterator(iter, ptr, size)
		//Iterator(OffsetsHolderIterator iter, PixelContainerPointer ptr ):ConstIterator(iter, ptr)
		{}

	private:
	};

	/** returns an iterator that points to the beginning of the container */
	Iterator Begin()
	{
		return Iterator( m_IdHolder.begin(), (PixelContainerPointer) m_Image->GetPixelContainer(), m_MeasurementVectorSize );
	}

	/** returns an iterator that points to the end of the container */
	Iterator End()
	{
		return Iterator( m_IdHolder.end(), (PixelContainerPointer) m_Image->GetPixelContainer(), m_MeasurementVectorSize  );
	}

	/** returns an iterator that points to the beginning of the container */
	ConstIterator Begin() const
	{
		return ConstIterator( m_IdHolder.begin(), (PixelContainerPointer) m_Image->GetPixelContainer(), m_MeasurementVectorSize  );
	}

	/** returns an iterator that points to the end of the container */
	ConstIterator End() const
	{
		return ConstIterator( m_IdHolder.end(), (PixelContainerPointer) m_Image->GetPixelContainer(), m_MeasurementVectorSize  );
	}

protected:
	MaskedImageToListSampleAdaptor();
	virtual ~MaskedImageToListSampleAdaptor() {}
	void PrintSelf(std::ostream & os, Indent indent) const;

private:
	MaskedImageToListSampleAdaptor(const Self &); //purposely not implemented
	void operator=(const Self &);//purposely not implemented

	ImageConstPointer m_Image;
	MaskConstPointer m_Mask;
	size_t m_MeasurementVectorSize;
	mutable TotalAbsoluteFrequencyType m_TotalFrequency;
	mutable MeasurementVectorType m_MeasurementVectorInternal;
	InstanceIdentifierHolder m_IdHolder;

}; // end of class MaskedImageToListSampleAdaptor
} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "MaskedImageToListSampleAdaptor.hxx"
#endif

#endif
