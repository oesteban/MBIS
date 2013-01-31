// --------------------------------------------------------------------------------------
// File:          ImageBiasedGMModelEMEstimator.h
// Date:          May 14, 2012
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

#ifndef IMAGEBIASEDGMMODELEMESTIMATOR_H_
#define IMAGEBIASEDGMMODELEMESTIMATOR_H_

#include <itkSimpleDataObjectDecorator.h>
#include <itkImage.h>
#include <itkMixtureModelComponentBase.h>
#include <itkLogVariableLengthVectorImageFilter.h>
#include <itkExpVariableLengthVectorImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>
#include "MaskedImageToListSampleAdaptor.h"
#include "GaussianMembershipFunction.h"
#include "BSplineRegressionImageFilter.h"
#include "GenerateModelImageFilter.h"
#include "NormalizeToModelImageFilter.h"

using namespace itk;

namespace mfbs {
namespace Statistics {

template<class TInputVectorImage, class TProbabilityPixelType = float >
class ImageBiasedGMModelEMEstimator: public itk::Object {
public:
	/** Standard class typedef */
	typedef ImageBiasedGMModelEMEstimator Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	/** Standard macros */
	itkTypeMacro(ImageBiasedGMModelEMEstimator, itk::Object);
	itkNewMacro (Self);


	typedef TInputVectorImage                                                      InputVectorImageType;
	typedef typename InputVectorImageType::Pointer                                 InputVectorPointer;
	typedef typename InputVectorImageType::ConstPointer                            InputVectorConstPointer;

    typedef TProbabilityPixelType                                                  ProbabilityPixelType;
    typedef itk::Image< ProbabilityPixelType,TInputVectorImage::ImageDimension >   ProbabilityImageType;
    typedef typename ProbabilityImageType::Pointer                                 ProbabilityImagePointer;
    typedef typename ProbabilityImageType::ConstPointer                            ProbabilityImageConstPointer;

    typedef typename std::vector< ProbabilityImagePointer >                        ProbabilityImageVector;

	typedef itk::Statistics::MaskedImageToListSampleAdaptor
	                    < InputVectorImageType, ProbabilityImageType >             InputSampleType;
	typedef typename InputSampleType::MeasurementType                              MeasurementType;
	typedef typename InputSampleType::MeasurementVectorType                        MeasurementVectorType;

/*
    typedef float                                                                  BiasPixelType;
    typedef itk::Image< BiasPixelType,
    		itk::GetImageDimension<TInputVectorImage>::ImageDimension >            BiasImageType;
    typedef typename BiasImageType::Pointer                                        BiasImagePointer;
    typedef typename BiasImageType::ConstPointer                                   BiasImageConstPointer;
    typedef typename std::vector< BiasImagePointer >                               BiasVector;

    typedef typename itk::MultiplyByConstantImageFilter
    		< ProbabilityImageType, MeasurementType, BiasImageType >               MultiplyFilter;
    		*/

    typedef typename itk::AddImageFilter
    		<InputVectorImageType, InputVectorImageType, InputVectorImageType>     AddFilter;
    typedef typename itk::SubtractImageFilter
    		<InputVectorImageType, InputVectorImageType, InputVectorImageType>     SubtractFilter;

    typedef typename mfbs::BSplineRegressionImageFilter<InputVectorImageType>      FieldEstimatorFilter;

    typedef typename mfbs::NormalizeToModelImageFilter<InputVectorImageType>       FieldNormalizerFilter;

    typedef typename mfbs::GenerateModelImageFilter< InputVectorImageType, ProbabilityImageType >
                                                                                   GenerateModelFilter;
    //typedef typename GenerateModelFilter::OutputImageType                        ModelImageType;
    //typedef typename ModelImageType::Pointer                                     ModelImagePointer;
    typedef typename GenerateModelFilter::InputChannelImageType                  ChannelImageType;

    typedef typename itk::LogImageFilter<InputVectorImageType, InputVectorImageType>         LogFilter;
    typedef typename itk::ExpImageFilter<InputVectorImageType, InputVectorImageType>         ExpFilter;




	/** Typedef requried to generate dataobject decorated output that can
	 * be plugged into SampleClassifierFilter */
	typedef GaussianMembershipFunction<MeasurementVectorType> GaussianMembershipFunctionType;
	typedef typename GaussianMembershipFunctionType::Pointer GaussianMembershipFunctionPointer;

	typedef itk::Statistics::MembershipFunctionBase<MeasurementVectorType> MembershipFunctionType;
	typedef typename MembershipFunctionType::ConstPointer MembershipFunctionPointer;
	typedef std::vector<MembershipFunctionPointer> MembershipFunctionVectorType;
	typedef itk::SimpleDataObjectDecorator<MembershipFunctionVectorType> MembershipFunctionVectorObjectType;
	typedef typename MembershipFunctionVectorObjectType::Pointer MembershipFunctionVectorObjectPointer;

	/** Type of the mixture model component base class */
	typedef itk::Statistics::MixtureModelComponentBase<InputSampleType> ComponentType;

	/** Type of the component pointer storage */
	typedef std::vector<ComponentType *> ComponentVectorType;

	/** Type of the membership function base class */
	typedef itk::Statistics::MembershipFunctionBase<MeasurementVectorType> ComponentMembershipFunctionType;

	/** Type of the array of the proportion values */
	typedef itk::Array<double> ProportionVectorType;

	itkSetConstObjectMacro( Input, InputVectorImageType );
	itkGetConstObjectMacro( Input, InputVectorImageType );

	itkSetConstObjectMacro( MaskImage, ProbabilityImageType );
	itkGetConstObjectMacro( MaskImage, ProbabilityImageType );

	/** Set/Gets the initial proportion values. The size of proportion
	 * vector should be same as the number of component (or classes) */
	void SetInitialProportions(ProportionVectorType & propotion);

	const ProportionVectorType & GetInitialProportions() const;

	/** Gets the result proportion values */
	const ProportionVectorType & GetProportions() const;

	/** typedef for decorated array of proportion */
	typedef itk::SimpleDataObjectDecorator<ProportionVectorType> MembershipFunctionsWeightsArrayObjectType;
	typedef typename MembershipFunctionsWeightsArrayObjectType::Pointer MembershipFunctionsWeightsArrayPointer;

	/** Get method for data decorated Membership functions weights array */
	const MembershipFunctionsWeightsArrayObjectType * GetMembershipFunctionsWeightsArray() const;

	/** Set/Gets the maximum number of iterations. When the optimization
	 * process reaches the maximum number of interations, even if the
	 * class parameters aren't converged, the optimization process
	 * stops. */
	void SetMaximumIteration(int numberOfIterations);

	int GetMaximumIteration() const;

	/** Gets the current iteration. */
	int GetCurrentIteration() {
		return m_CurrentIteration;
	}

	/** Adds a new component (or class). */
	int AddComponent(ComponentType *component);

	/** Gets the total number of classes currently plugged in. */
	unsigned int GetNumberOfComponents() const;

	/** Runs the optimization process. */
	void Update();

	/** Termination status after running optimization */
	enum TERMINATION_CODE {
		CONVERGED = 0, NOT_CONVERGED = 1
	};

	/** Gets the termination status */
	TERMINATION_CODE GetTerminationCode() const;

	/** Gets the membership function specified by componentIndex
	 argument. */
	ComponentMembershipFunctionType * GetComponentMembershipFunction(int componentIndex) const;

	const ProbabilityImageType* GetPosterior( const size_t id ) const {
		return m_Posteriors[id];
	}

	/** Output Membership function vector containing the membership functions with
	 * the final optimized parameters */
	const MembershipFunctionVectorObjectType * GetOutput() const;

	itkGetConstObjectMacro( CorrectedInput, InputVectorImageType );
	itkGetConstObjectMacro( CurrentBias, InputVectorImageType );

	itkGetConstMacro(UseBiasCorrection, bool);
	itkSetMacro(UseBiasCorrection, bool);

protected:
	ImageBiasedGMModelEMEstimator();
	virtual ~ImageBiasedGMModelEMEstimator() {
	}
	void PrintSelf(std::ostream & os, Indent indent) const;

	bool CalculateDensities();

	double CalculateExpectation() const;

	bool UpdateComponentParameters();

	bool UpdateProportions();

	bool UpdateBiasFieldEstimate();

	/** Starts the estimation process */
	void GenerateData();

private:
	/** Target data sample pointer*/
	typename InputSampleType::Pointer         m_Sample;

	int m_MaxIteration;
	int m_CurrentIteration;

	TERMINATION_CODE m_TerminationCode;
	ComponentVectorType m_ComponentVector;
	ProportionVectorType m_InitialProportions;
	ProportionVectorType m_Proportions;
	ProbabilityImageVector m_Posteriors;

	//std::vector<std::vector<double> > m_Ptly;

	MembershipFunctionVectorObjectPointer m_MembershipFunctionsObject;
	MembershipFunctionsWeightsArrayPointer m_MembershipFunctionsWeightArrayObject;


	InputVectorConstPointer      m_Input;
	InputVectorPointer           m_CorrectedInput;
	ProbabilityImageConstPointer m_MaskImage;
	InputVectorPointer           m_BiasLog;
	InputVectorPointer           m_CurrentBias;
	bool                         m_UseBiasCorrection;
	bool                         m_BiasCorrectionStopped;

}; // end of class

}// end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "ImageBiasedGMModelEMEstimator.hxx"
#endif


#endif /* IMAGEBIASEDGMMODELEMESTIMATOR_H_ */
