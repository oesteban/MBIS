// --------------------------------------------------------------------------------------
// File:    MultispectralMRFFilter.h
// Date:    Nov 9, 2011
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

#ifndef MULTISPECTRALMRFFILTER_H_
#define MULTISPECTRALMRFFILTER_H_

#include <vector>

#include <itkConstNeighborhoodIterator.h>
#include <itkNeighborhoodIterator.h>
#include <itkTimeProbesCollectorBase.h>

#include "MultispectralFilter.h"
//#include "MultispectralGaussianEstimator.h"
#include "MultichannelParameterEstimator.h"

#include "MRFEnergy.h"
#include "MRFOptimizer.h"
#include "MRFExpansionOptimizer.h"
#include "MRFSwapOptimizer.h"




namespace itk
{

template <class TInputComponent, class TProbabilityPixelType = float, class TEnergyValType = float >
class ITK_EXPORT MultispectralMRFFilter:
	public MultispectralFilter< TInputComponent, itk::Image< unsigned char, TInputComponent::ImageDimension> >
{
public:
	/** Standard class typedefs */
	typedef MultispectralMRFFilter                                                 Self;
	typedef SmartPointer<Self>                                                     Pointer;
	typedef SmartPointer<const Self>                                               ConstPointer;

	typedef TProbabilityPixelType                                                  RealType;
	typedef TEnergyValType                                                         EnergyValType;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);
	/** Run-time type information (and related methods). */
	itkTypeMacro(MultispectralMRFFilter, MultispectralFilter);

	/** Dimension of the images. */
	itkStaticConstMacro( ImageDimension, size_t, TInputComponent::ImageDimension );

	typedef unsigned char                                                          LabelType;
	typedef Image< LabelType, ImageDimension >                                     OutputImageType;
	typedef typename OutputImageType::Pointer                                      OutputImagePointer;
	typedef typename itk::ImageRegionIterator< OutputImageType >                   OutputImageIterator;
	typedef typename itk::ImageRegionConstIterator< OutputImageType >              OutputImageConstIterator;
	typedef typename itk::ConstNeighborhoodIterator<OutputImageType>               ConstLabelNeighIterator;

	typedef MultispectralFilter< TInputComponent, OutputImageType >                Superclass;

	typedef typename Superclass::InputVectorImageType                              InputVectorImageType;
	typedef typename InputVectorImageType::ConstPointer                            InputVectorConstPointer;
	typedef typename InputVectorImageType::PixelType		      				   MeasurementVectorType;

    typedef typename Superclass::ParametersType                                    ParametersType;
    typedef typename Superclass::ParametersVectorType                              ParametersVectorType;

    /** Priors probability images typedefs */
    typedef TProbabilityPixelType                                                  ProbabilityPixelType;
    typedef itk::Image< ProbabilityPixelType, ImageDimension >                     ProbabilityImageType;
    typedef typename ProbabilityImageType::Pointer                                 ProbabilityImagePointer;
    typedef typename ProbabilityImageType::ConstPointer                            ProbabilityImageConstPointer;
    typedef typename itk::ImageRegionConstIterator< ProbabilityImageType >         ProbabilityImageConstIterator;
    typedef typename std::vector< ProbabilityImageConstPointer >                   ProbabilityConstImagesVector;
    typedef typename std::vector< ProbabilityImagePointer >                        ProbabilityImagesVector;
    typedef typename itk::ConstNeighborhoodIterator<ProbabilityImageType>          ConstProbNeighIterator;
    //typedef itk::MultispectralGaussianEstimator< ProbabilityImageType >            GaussianEstimatorFilter;

    typedef itk::VariableSizeMatrix< EnergyValType >                               EnergyMatrixType;
    typedef itk::Image< EnergyValType, ImageDimension >                            EnergyMapType;

	typedef itk::MultichannelParameterEstimator< TInputComponent >                 ParameterEstimator;



	typedef itk::Statistics::MaskedImageToListSampleAdaptor
	                    < InputVectorImageType, ProbabilityImageType >             InputSampleType;

    typedef typename  itk::Statistics::WeightedCovarianceSampleFilter
    		                                    < InputSampleType >                CovarianceEstimatorType;
    typedef typename CovarianceEstimatorType::WeightArrayType                      WeightArrayType;

    typedef typename mfbs::Statistics::GaussianMixtureModelComponent
    		                                     < InputSampleType >               ComponentType;
    typedef std::vector< typename ComponentType::Pointer >                         ComponentVectorType;

    typedef FixedArray<unsigned, ImageDimension>        ArrayType;

    typedef MRFEnergy
    	< InputVectorImageType, ProbabilityImageType, LabelType, EnergyValType >   EnergyType;

	typedef MRFOptimizer< EnergyType >                                             OptimizerType;
	typedef MRFGCOptimizer< EnergyType >                                           GCOptimizerType;
	typedef SwapOptimizer< EnergyType >                                            SwapOptimizerType;
	typedef ExpansionOptimizer< EnergyType >                                       ExpansionOptimizerType;


    itkGetConstMacro( NumberOfClasses, unsigned char);

    void SetNumberOfClasses( unsigned char c ) {
    	m_NumberOfClasses = c;
    }

    void SetInitialParameters( ParametersVectorType& params) { m_InitialParameters = params; }
    itkGetConstMacro( InitialParameters, ParametersVectorType );

    itkSetObjectMacro( MaskImage, ProbabilityImageType);
    itkGetConstObjectMacro( MaskImage, ProbabilityImageType);

    itkSetObjectMacro( InitialLabeling, OutputImageType );
    itkGetConstObjectMacro( InitialLabeling, OutputImageType );

    itkGetConstObjectMacro( Input, InputVectorImageType);

    itkGetConstMacro( EnergyMatrix, EnergyMatrixType );
    itkGetConstMacro( MatrixOffset, size_t );

    itkSetMacro( Lambda, double );
    itkGetConstMacro( Lambda, double );

    itkSetMacro( ExternalLambda, double );
    itkGetConstMacro( ExternalLambda, double );

    itkSetClampMacro( MRFIterations, size_t, 1, 100 );
    itkGetConstMacro( MRFIterations, size_t );

    itkSetMacro(UseExplicitPVModel, bool);
    itkGetConstMacro(UseExplicitPVModel, bool);


    /** Minimization algorithm to be used. */
    virtual void MinimizeFunctional();

    void ExtendToExplicitPVModel();

    void WriteIteration( );
    void WriteMaps( );

    const ProbabilityImageType* GetPosteriorProbabilityImage( size_t label ) {
    	return m_EnergyFunction->GetPostProbMapForLabel(label);
    }

    itkSetMacro( UseOutputPriors, bool );
    itkGetConstMacro( UseOutputPriors, bool );

    itkSetMacro( UseCachedEnergyMap, bool );
    itkGetConstMacro( UseCachedEnergyMap, bool );

    itkSetMacro( MinimizationMode, size_t );
    itkGetConstMacro( MinimizationMode, size_t );

    itkSetMacro( OutputPrefix, std::string );
    itkGetConstMacro( OutputPrefix, std::string );
/*
    ComponentType *  GetComponent( size_t i ) { return m_CVector[i]; }
	ComponentVectorType& GetComponentVector() { return m_CVector; }

    void SetComponent( size_t i, ComponentType* comp ) { m_CVector[i] = comp; }
	void SetComponentVector( ComponentVectorType& v ) { m_CVector = v; }*/

	itkSetMacro( MatrixFile , std::string );
	itkGetConstMacro( MatrixFile, std::string );

	void  SetSpatialEnergyMap( LabelType l, EnergyMapType* m ) {
		m_SpatialEnergyMapIndex.push_back(l);
		m_SpatialEnergyMap.push_back(m);
	}


protected:
	MultispectralMRFFilter();
	virtual ~MultispectralMRFFilter();

	virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;


	void InitializeModel();
	void InitializeExternalFieldMaps();
	void LoadMatrixFile();

	/** Standard pipeline method. While this class does not implement a
	 * ThreadedGenerateData(), its GenerateData() delegates all
	 * calculations to an InputSampleFilter. Multithreading depends on
	 * that one  */
	void GenerateData();

private:
	MultispectralMRFFilter( const Self& ); //purposely not implemented
	void operator=(const Self& );       //purposely not implemented

	InputVectorConstPointer                                    m_Input;
	typename InputSampleType::Pointer                          m_Sample;
	unsigned char                                              m_NumberOfComponents;
	unsigned char                                              m_NumberOfClasses;
	size_t                                                     m_MatrixOffset;
	unsigned char                                              m_NumberOfPureTissues;
	size_t                                                     m_MRFIterations;
	size_t													   m_CurrentIteration;

	ComponentVectorType                                        m_CVector;
	typename EnergyType::Pointer                               m_EnergyFunction;

	ProbabilityImageConstPointer                               m_MaskImage;
	std::vector< LabelType >                                   m_SpatialEnergyMapIndex;
	ProbabilityImagesVector                                    m_SpatialEnergyMap;
	ParametersVectorType                                       m_InitialParameters;
	RealType                                                   m_Lambda;
	RealType                                                   m_ExternalLambda;
	EnergyMatrixType                                           m_EnergyMatrix;
	bool                                                       m_UseOutputPriors;
	bool                                                       m_UseExplicitPVModel;
	bool                                                       m_UseCachedEnergyMap;


	std::vector< ProbabilityImagePointer >                     m_GaussianProbabilityMap;

	typename NeighborhoodIterator<OutputImageType>::RadiusType m_Radius;


	//typename GaussianEstimatorFilter::Pointer                  m_GaussianEstimator;
	OutputImagePointer                                         m_InitialLabeling;
	size_t                                                     m_MinimizationMode;

	std::string                                                m_MatrixFile;
	std::string                                                m_OutputPrefix;


};

}


#include "MultispectralMRFFilter.txx"


#endif /* MULTISPECTRALMRFFILTER_H_ */
