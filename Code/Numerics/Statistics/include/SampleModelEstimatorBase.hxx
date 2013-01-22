// --------------------------------------------------------------------------------------
// File:          SampleModelEstimatorBase.hxx
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

#ifndef SAMPLEMODELESTIMATORBASE_HXX_
#define SAMPLEMODELESTIMATORBASE_HXX_

namespace itk {

template< class TInputSample, class TOutputType >
SampleModelEstimatorBase<TInputSample,TOutputType>
::SampleModelEstimatorBase() {
	// Create the output
	typename TOutputType::Pointer output = static_cast< TOutputType* >( this->MakeOutput(0).GetPointer() );
	this->ProcessObject::SetNumberOfRequiredOutputs(1);
	this->ProcessObject::SetNthOutput( 0, output.GetPointer() );

	// Set default behavior
	this->ReleaseDataBeforeUpdateFlagOff();
}

template< class TInputSample, class TOutputType >
ProcessObject::DataObjectPointer
SampleModelEstimatorBase<TInputSample,TOutputType>
::MakeOutput( ProcessObject::DataObjectPointerArraySizeType ){
	return TOutputType::New().GetPointer();
}


template< class TInputSample, class TOutputType >
typename SampleModelEstimatorBase<TInputSample,TOutputType>::OutputType *
SampleModelEstimatorBase<TInputSample,TOutputType>
::GetOutput() {
	return itkDynamicCastInDebugMode< TOutputType * >( this->GetPrimaryOutput() );
}

template< class TInputSample, class TOutputType >
const typename SampleModelEstimatorBase<TInputSample,TOutputType>::OutputType *
SampleModelEstimatorBase<TInputSample,TOutputType>
::GetOutput() const {
	return itkDynamicCastInDebugMode< const TOutputType * >( this->GetPrimaryOutput() );
}

template< class TInputSample, class TOutputType >
typename SampleModelEstimatorBase<TInputSample,TOutputType>::OutputType *
SampleModelEstimatorBase<TInputSample,TOutputType>
::GetOutput(size_t idx) {
	TOutputType *out = dynamic_cast< TOutputType* >( this->ProcessObject::GetOutput(idx) );

	if ( out == NULL & this->ProcessObject::GetOutput(idx) != NULL ) {
		itkWarningMacro( << "Unable to convert output number " << idx << " to type " << typeid( OutputType ).name () );
	}
	return out;
}

template< class TInputSample, class TOutputType >
void
SampleModelEstimatorBase<TInputSample,TOutputType>
::GraftOutput(const DataObjectIdentifierType & key, DataObject *graft) {
	if (! graft) {
		itkExceptionMacro(<< "Requested to graft output that is a NULL pointer");
	}

	// we use the process object method since all out output may not be
	// of the same type
	DataObject *output = this->ProcessObject::GetOutput(key);

	// Call GraftImage to copy meta-information, regions, and the pixel container
	output->Graft(graft);
}


template< class TInputSample, class TOutputType >
void
SampleModelEstimatorBase<TInputSample,TOutputType>
::GraftNthOutput(size_t idx, DataObject *graft ) {
	if ( idx >= this->GetNumberOfIndexedOutputs() ) {
		itkExceptionMacro(<< "Requested to graft output " << idx
				<< " but this filter only has " << this->GetNumberOfIndexedOutputs() << " indexed Outputs.");
	}
	this->GraftOutput( this->MakeNameFromOutputIndex(idx), graft );
}

//-------------------------------------------------------------------------------------
template< class TInputSample, class TOutputType >
void
SampleModelEstimatorBase<TInputSample,TOutputType>
::GenerateData() {
	this->AllocateOutputs();
	this->BeforeThreadedGenerateData();

	ThreadStruct str;
	str.Filter = this;

	this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
	this->GetMultiThreader()->SetSingleMethod( this->ThreaderCallback, &str );
	this->GetMultiThreader()->SingleMethodExecute();
	this->AfterThreadedGenerateData();
}

template< class TInputSample, class TOutputType >
void
SampleModelEstimatorBase<TInputSample,TOutputType>
::ThreadedGenerateData(const InputSubSampleType* subSample, ThreadIdType threadId ){
	// The following code is equivalent to:
	// itkExceptionMacro("subclass should override this method!!!");
	// The ExceptionMacro is not used because gcc warns that a
	// 'noreturn' function does return
	std::ostringstream message;

	message << "itk::ERROR: " << this->GetNameOfClass()
	        << "(" << this << "): " << "Subclass should override this method!!!" << std::endl
	        << "The signature of ThreadedGenerateData() has been changed in ITK v4 to use the new ThreadIdType." << std::endl
	        << this->GetNameOfClass() << "::ThreadedGenerateData() might need to be updated to used it.";
	ExceptionObject e_(__FILE__, __LINE__, message.str().c_str(), ITK_LOCATION);
	throw e_;
}

template< class TInputSample, class TOutputType >
ITK_THREAD_RETURN_TYPE
SampleModelEstimatorBase<TInputSample,TOutputType>
::ThreaderCallback( void *arg ) {
	ThreadStruct *str;
	ThreadIdType total, threadId, threadCount;

	threadId = ((MultiThreader::ThreadInfoStruct *) (arg))->ThreadID;
	threadCount = ((MultiThreader::ThreadInfoStruct *) (arg))->NumberOfThreads;

	str = (ThreadStruct *) (((MultiThreader::ThreadInfoStruct *) (arg))->UserData);

	// execute the actual method with appropriate output region
	// first find out how many pieces extent can be split into.
	InputSubSamplePointer subSample;
	total = str->Filter->SplitRequestedRegion(threadId, threadCount, subSample);

	if (threadId < total) {
		str->Filter->ThreadedGenerateData(subSample, threadId);
	}
	// else
	//   {
	//   otherwise don't use this thread. Sometimes the threads dont
	//   break up very well and it is just as efficient to leave a
	//   few threads idle.
	//   }

	return ITK_THREAD_RETURN_VALUE;
}


//----------------------------------------------------------------------------
template< class TInputSample, class TOutputType >
size_t
SampleModelEstimatorBase<TInputSample,TOutputType>
::SplitRequestedRegion(size_t i, size_t num, InputSubSampleType & subSample) {
	// Split the region


	// set the split region ivars
	subSample.SetIndex(splitIndex);
	subSample.SetSize(splitSize);

	itkDebugMacro("  Split Piece: " << subSample);

	return maxThreadIdUsed + 1;
}

} // end namespace itk

#endif /* SAMPLEMODELESTIMATORBASE_HXX_ */
