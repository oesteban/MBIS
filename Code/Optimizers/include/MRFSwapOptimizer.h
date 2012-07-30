// --------------------------------------------------------------------------------------
// File:    MRFSwapOptimizer.h
// Date:    Mar 15, 2012
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

#ifndef MRFSWAPOPTIMIZER_H_
#define MRFSWAPOPTIMIZER_H_

#include <vector>
#include "MRFGCOptimizer.h"

template < class TEnergy >
class SwapOptimizer:
	public MRFGCOptimizer< TEnergy >
	{
public:
	typedef SwapOptimizer< TEnergy >                Self;
	typedef itk::SmartPointer<Self>                 Pointer;
	typedef itk::SmartPointer<const Self>           ConstPointer;
	typedef MRFGCOptimizer< TEnergy >               Super;


	typedef typename Super::EnergyValType            EnergyValType;
	typedef typename Super::EnergyType               EnergyType;
	typedef typename Super::VarVector                VarVector;
	typedef typename Super::SiteVector               SiteVector;
	typedef typename Super::LookupTable              LookupTable;
	typedef typename Super::LabelType                LabelType;
	typedef typename Super::GraphType                GraphType;

	itkStaticConstMacro( ImageDimension, size_t, TEnergy::ImageDimension );
	typedef itk::Image< LabelType, ImageDimension >                        LabelImageType;
	typedef typename LabelImageType::Pointer                               LabelImagePointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);
	/** Run-time type information (and related methods). */
	itkTypeMacro(SwapOptimizer, MRFGCOptimizer);

protected:
	SwapOptimizer();
	virtual ~SwapOptimizer() {};


	void PerformMove() {
		// 3. For each pair of labels {alpha,beta} do swap
		for (size_t l1Id = 0;  l1Id < this->m_NumberOfClasses;  l1Id++ ) {
			for (int l2Id = this->m_NumberOfClasses-1;  l2Id >= 0 ;  l2Id-- ) {
				size_t label1 = this->m_LabelTable[l1Id];
				size_t label2 = this->m_LabelTable[l2Id];

				if( label1 < label2 )	PerformSwap(label1, label2);
			}
		}
	}
private:
	SwapOptimizer( const Self& ); //purposely not implemented
	void operator=(const Self& );       //purposely not implemented


	void PerformSwap( LabelType aLabel, LabelType bLabel );
	void SetUpDataCostsSwap( LabelType aLabel, LabelType bLabel, GraphType *e );
	void SetUpSmoothCostsSwap( LabelType aLabel, LabelType bLabel, GraphType *e );

};


#include "MRFSwapOptimizer.txx"

#endif /* MRFSWAPOPTIMIZER_H_ */
