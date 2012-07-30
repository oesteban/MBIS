// --------------------------------------------------------------------------------------
// File:          MRFGCOptimizer.h
// Date:          Mar 15, 2012
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

#ifndef MRFGCOPTIMIZER_H_
#define MRFGCOPTIMIZER_H_

#ifndef GCO_MAX_ENERGYTERM          // maximum safe coefficient to avoid integer overflow
#define GCO_MAX_ENERGYTERM 1000000  // if a data/smooth/label cost term is larger than this,
#endif                              // the library will raise an exception


#include <energy.h>
#include "MRFOptimizer.h"

template < class TEnergy >
class MRFGCOptimizer: public MRFOptimizer<TEnergy> {

public:

	typedef MRFGCOptimizer<TEnergy>                                        Self;
	typedef itk::SmartPointer<Self>                                        Pointer;
	typedef itk::SmartPointer<const Self>                                  ConstPointer;

	typedef MRFOptimizer<TEnergy>                                          Super;
	typedef typename Super::EnergyType                                     EnergyType;
	typedef typename EnergyType::SiteVector                                SiteVector;
	typedef typename EnergyType::EnergyValType                             EnergyValType;
	typedef typename EnergyType::LabelType                                 LabelType;
	typedef typename EnergyType::LabelImageType                            LabelImageType;
	typedef itk::Array< LabelType >                                        LabelTableType;

	typedef Energy<EnergyValType, EnergyValType, EnergyValType>            GraphType;
	typedef std::vector< typename GraphType::Var >                         VarVector;

	typedef std::vector< int >                                             LookupTable;

	enum MINIMIZATION_MODES {
		EXPANSION,
		SWAP
	};

	/** Run-time type information (and related methods). */
	itkTypeMacro(MRFGCOptimizer, MRFOptimizer);

    void StartOptimization();

    itkSetClampMacro( MaxAlgorithmIterations, size_t, 1, 100 );
    itkGetConstMacro( MaxAlgorithmIterations, size_t );
    itkGetConstMacro( LabelTable, LabelTableType );

protected:
    MRFGCOptimizer();
	~MRFGCOptimizer() {};

    void AddtLink( GraphType *e, size_t i, EnergyValType e0, EnergyValType e1);
    void AddnLink( GraphType *e, size_t i, size_t j, EnergyValType e00, EnergyValType e01, EnergyValType e10, EnergyValType e11 );
    void ScrambleLabelTable();
    inline void ResetLookupSiteVar() {
		std::fill(this->m_LookupSiteVar.begin(),this->m_LookupSiteVar.end(), -1 );
	}
    void optimize();

    //virtual inline bool activeCondition( LabelType label ) const = 0;

    virtual void PerformMove() = 0;

	typename LabelImageType::Pointer  m_Output;
	EnergyValType                     m_BeforeMoveEnergy;
	SiteVector                        m_ActiveSites;
	LookupTable                       m_LookupSiteVar;
	LabelTableType                    m_LabelTable;
	size_t                            m_NumberOfClasses;

private:
	size_t                            m_MaxAlgorithmIterations;

};

#include "MRFGCOptimizer.txx"


#endif /* MRFGCOPTIMIZER_H_ */
