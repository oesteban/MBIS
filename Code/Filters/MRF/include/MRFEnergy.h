// --------------------------------------------------------------------------------------
// File:    MRFEnergy.h
// Date:    Mar 16, 2012
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

#ifndef MRFENERGY_H_
#define MRFENERGY_H_


#ifndef MRF_INF                     // maximum safe coefficient to avoid integer overflow
#define MRF_INF 100000              // if a data/smooth/label cost term is larger than this,
#endif                              // the library will raise an exception

#ifndef DATA_MAX_ENERGY
#define DATA_MAX_ENERGY 100
#endif

#include <vector>

#include <itkObject.h>
#include <itkImage.h>
#include "Numerics/Statistics/include/MaskedImageToListSampleAdaptor.h"
#include "Numerics/Statistics/include/MembershipEnergyFunctionBase.h"
#include "Numerics/Statistics/include/GaussianMembershipFunction.h"

template < class TInputImage, class TMaskImage, class TLabelType = unsigned char, class TEnergyVal = float >
class MRFEnergy : public itk::Object {
public:
	typedef MRFEnergy<TInputImage, TMaskImage, TLabelType, TEnergyVal>     Self;
	typedef itk::SmartPointer<Self>                                        Pointer;
	typedef itk::SmartPointer<const Self>                                  ConstPointer;

	typedef TInputImage                                                    InputImageType;
	typedef typename InputImageType::ConstPointer                          InputImageConstPointer;
	typedef typename InputImageType::PixelType		      			       MeasurementVectorType;

	itkStaticConstMacro( ImageDimension, size_t, TInputImage::ImageDimension );

	typedef TMaskImage                                                     MaskImageType;
	typedef typename MaskImageType::ConstPointer                           MaskImageConstPointer;

	typedef itk::Statistics::MaskedImageToListSampleAdaptor
	                    < InputImageType, MaskImageType >                  InputSampleType;

	typedef TLabelType                                                     LabelType;
	typedef itk::Image< LabelType, ImageDimension >                        LabelImageType;
	typedef typename LabelImageType::Pointer                               LabelImagePointer;

	typedef TEnergyVal                                                     EnergyValType;
	typedef itk::Image< EnergyValType, ImageDimension >                    EnergyMapType;
	typedef std::vector< typename EnergyMapType::Pointer >                 EnergyMapVector;
	typedef itk::VariableSizeMatrix< EnergyValType >                       EnergyMatrixType;

	typedef std::vector< size_t >                                          SiteVector;


    typedef typename mfbs::Statistics
    		::MembershipEnergyFunctionBase< MeasurementVectorType >        MembershipFunction;
    typedef std::vector< typename MembershipFunction::Pointer >            MembershipFunctionVectorType;
    typedef typename mfbs::Statistics
        		::GaussianMembershipFunction< MeasurementVectorType >      DefaultMembershipFunction;

    typedef itk::VariableSizeMatrix< double >                              MatrixType;
    typedef itk::Array< double >                                           ParametersType;


	/** Method for creation through the object factory. */
	itkNewMacro(Self);
	/** Run-time type information (and related methods). */
	itkTypeMacro(MRFEnergy, Object);


	/** Pixel-wise energy calculators **/
    inline EnergyValType GetSampleEnergyForLabel( LabelType label, size_t offset) {
		assert( label > 0 );
    	EnergyValType e = *(this->m_DataEnergyMapBuffer[label-1] + offset) +
    			*(this->m_SpatialEnergyMapBuffer[label-1] + offset);

    	if (e>MRF_INF) return MRF_INF;
   		return e;
    }
    inline EnergyValType GetSampleEnergy( size_t offset) {
    	LabelType label = *(this->m_OutputBuffer + offset );
   		return this->GetSampleEnergyForLabel( label, offset );
    }
    inline EnergyValType GetSampleEnergyWithIndex( LabelType label, typename LabelImageType::IndexType &siteIdx) {
   		return this->m_DataEnergyMap[label-1]->GetPixel( siteIdx ) + this->m_SpatialEnergyMap[label-1]->GetPixel( siteIdx );
    }
	inline EnergyValType GetNeighborCost( LabelType l1, LabelType l2 ) {
    	size_t increment = (this->m_MatrixOffset==0)?1:0;
    	EnergyValType val = this->m_EnergyMatrix(l1 - increment,l2 - increment);

    	return (val==MRF_INF)?MRF_INF:val*m_Lambda;
    }

	EnergyValType GetCliqueEnergy(size_t offset);
	EnergyValType GetCliqueEnergyForLabel(size_t offset, LabelType cLabel);
	EnergyValType GetExternalFieldEnergy(size_t offset);
	EnergyValType GetExternalFieldEnergyForLabel(size_t offset, LabelType cLabel);



    itkGetObjectMacro( Output, LabelImageType );
    itkSetObjectMacro( Output, LabelImageType );

    itkGetConstObjectMacro( Input, InputImageType );
    itkSetConstObjectMacro( Input, InputImageType );

    itkGetConstObjectMacro( MaskImage, MaskImageType );
    itkSetConstObjectMacro( MaskImage, MaskImageType );

    itkSetClampMacro( Lambda, EnergyValType, 0.0, 1000.0 );
    itkGetConstMacro( Lambda, EnergyValType );

    itkSetClampMacro( ExternalLambda, EnergyValType, 0.0, 1000.0 );
    itkGetConstMacro( ExternalLambda, EnergyValType );

    itkSetClampMacro( NumberOfClasses, size_t, 2, 1000 );
    itkGetConstMacro( NumberOfClasses, size_t );

    itkGetConstMacro( EnergyMatrix, EnergyMatrixType );

    void SetEnergyMatrix( EnergyMatrixType m ) {
    	for ( size_t i = 0; i < m.Cols(); i++ ) {
    		for( size_t j = 0; j < m.Rows(); j++ ) {
    			if( m(i,j) < MRF_INF ) m(i,j) * m_Lambda;
    			else
    				m(i,j) = MRF_INF;
    		}
    	}
    	m_EnergyMatrix = m;
    }

    const SiteVector GetIdHolder() const { return m_IdHolder; }

    LabelType* GetOutputPointer() { return m_Output->GetBufferPointer(); }

    inline EnergyValType GetTotalDataEnergy() {
    	if( m_ModelModified ) UpdateMaps();
    	else if ( m_LabelingModified ) UpdateTotalDataEnergy();
    	return this->m_TotalDataEnergy;
    }

    inline EnergyValType GetTotalSmoothnessEnergy() {
    	if( m_LabelingModified ) UpdateTotalSmoothnessEnergy();
    	return this->m_TotalSmoothnessEnergy;
    }

    inline  EnergyValType GetTotalEnergy() { return GetTotalDataEnergy() + GetTotalSmoothnessEnergy(); }


    const EnergyMapType* GetDataEnergyMapForLabel( LabelType l )    { return m_DataEnergyMap[l-1]; }
    const EnergyMapType* GetDataProbMapForLabel  ( LabelType l )    { return m_DataProbMap[l-1]; }
    const EnergyMapType* GetPostProbMapForLabel  ( LabelType l )    { return m_PostProbMap[l-1]; }
    const EnergyMapType* GetSpatialEnergyMapForLabel( LabelType l ) { return m_SpatialEnergyMap[l-1]; }
    const EnergyMapType* GetMRFProbMapForLabel   ( LabelType l )    { return m_MRFProbMap[l-1]; }

    void SetSpatialEnergyMap( LabelType l, EnergyMapType* m );
    void SetMembershipFunction( size_t i, MembershipFunction* comp ) { m_MembershipVector[i] = comp; }
    void SetModelParametersForLabel( LabelType l, ParametersType p );

	void Initialize();
    void UpdateMaps();
    void ComputePostProbMaps();

    inline void ModelModified() { m_ModelModified = true; }
    inline void LabelingModified() { m_LabelingModified = true; }

    void SaveMRFEnergyMap( std::string name );

protected:
	MRFEnergy();
	virtual ~MRFEnergy() {};

private:

    /* returns the data part of the energy */
    void UpdateTotalDataEnergy();

    /* returns the smoothness part of the energy */
    void UpdateTotalSmoothnessEnergy();


	size_t                             m_NumberOfClasses;
	EnergyValType                      m_Lambda;
	EnergyValType                      m_ExternalLambda;

	InputImageConstPointer             m_Input;
	MaskImageConstPointer              m_MaskImage;
	LabelImagePointer                  m_Output;
	typename InputSampleType::Pointer  m_Sample;


	EnergyMatrixType                   m_EnergyMatrix;
	size_t                             m_MatrixOffset;

	bool                               m_IsInitialized;

	MembershipFunctionVectorType       m_MembershipVector;

	EnergyMapVector	                   m_DataEnergyMap;
	EnergyMapVector	                   m_DataProbMap;
	EnergyMapVector                    m_PostProbMap;
	EnergyMapVector                    m_MRFProbMap;
	std::vector < typename EnergyMapType::Pointer >  m_SpatialEnergyMap;

	std::vector < EnergyValType* >     m_DataEnergyMapBuffer;
	std::vector < EnergyValType* >     m_DataProbMapBuffer;
	std::vector < EnergyValType* >     m_PostProbMapBuffer;
	std::vector < EnergyValType* >     m_SpatialEnergyMapBuffer;
	std::vector < EnergyValType* >     m_MRFProbMapBuffer;

	EnergyValType                      m_TotalDataEnergy;
	EnergyValType                      m_TotalSmoothnessEnergy;
	LabelType*                         m_OutputBuffer;
	SiteVector                         m_IdHolder;

	bool m_ModelModified;
	bool m_LabelingModified;

	static const float MAX_LOG;
};

#include "MRFEnergy.txx"

#endif /* MRFENERGY_H_ */
