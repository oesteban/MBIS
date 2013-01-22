// --------------------------------------------------------------------------------------
// File:          SampleModelEstimatorBase.h
// Date:          Jan 22, 2013
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

#ifndef SAMPLEMODELESTIMATORBASE_H_
#define SAMPLEMODELESTIMATORBASE_H_

#include <itkProcessObject.h>
#include <itkSample.h>
#include <itkSubsample.h>

namespace itk {

/** \class SampleModelEstimatorBase
 *  \brief SampleModelEstimatorBase is the base class for all process objects that output
 *  some form of parameters and require a sample as input. Specifically, this class defines
 *  the SetInput() method for defining the input to a filter.
 *
 *  This class provides the infrastructure for supporting multithreaded processing of samples.
 *
 *  \ingroup SampleFilters
 *  \ingroup ITKCommon
 */

template< class TInputSample, class TOutputType >
class ITK_EXPORT SampleModelEstimatorBase: public ProcessObject {
public:
	typedef SampleModelEstimatorBase                     Self;
	typedef ProcessObject                                Superclass;
	typedef SmartPointer< Self >                         Pointer;
	typedef SmartPointer< const Self >                   ConstPointer;


	/** Smart Pointer type to a DataObject. */
	typedef DataObject::Pointer                          DataObjectPointer;
	typedef Superclass::DataObjectIdentifierType         DataObjectIdentifierType;
	typedef Superclass::DataObjectPointerArraySizeType   DataObjectPointerArraySizeType;

	/** Run-time type information (and related methods). */
	itkTypeMacro(SampleModelEstimatorBase, ProcessObject);

	typedef TInputSample                                 InputSampleType;
	typedef typename InputSampleType::Pointer            InputSamplePointer;
	typedef typename itk::Statistics::
			             Subsample<InputSampleType>      InputSubSampleType;
	typedef typename InputSubSampleType::Pointer         InputSubSamplePointer;

	typedef TOutputType                                  OutputType;


	OutputType * GetOutput(void);
	const OutputType * GetOutput(void) const;
	OutputType * GetOutput(size_t idx);

	virtual void GraftOutput( DataObject *output );
	virtual void GraftOutput( const DataObjectIdentifierType & key, DataObject *output );
	virtual void GraftNthOutput(size_t idx, DataObject *output);

	using Superclass::MakeOutput;
	virtual ProcessObject::DataObjectPointer MakeOutput(ProcessObject::DataObjectPointerArraySizeType idx);

protected:
	SampleModelEstimatorBase();
	virtual ~SampleModelEstimatorBase() {}

	virtual void GenerateData();
	virtual void ThreadedGenerateData(const InputSubSampleType* subSample, ThreadIdType threadId );
	virtual void BeforeThreadedGenerateData() {}
	virtual void AfterThreadedGenerateData() {}

	static ITK_THREAD_RETURN_TYPE ThreaderCallback(void *arg);

	struct ThreadStruct {
		Pointer Filter;
	};

private:
	SampleModelEstimatorBase( const Self &); // purposely not implemented
	void operator=(const Self &);            // purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "SampleModelEstimatorBase.hxx"
#endif


#endif /* SAMPLEMODELESTIMATORBASE_H_ */
