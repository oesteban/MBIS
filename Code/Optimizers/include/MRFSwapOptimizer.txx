// --------------------------------------------------------------------------------------
// File:          MRFSwapOptimizer.txx
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

#ifndef MRFSWAPOPTIMIZER_TXX_
#define MRFSWAPOPTIMIZER_TXX_

#include "MRFSwapOptimizer.h"

template < class TEnergy >
SwapOptimizer< TEnergy >::
SwapOptimizer() {

};

template < class TEnergy >
void
SwapOptimizer< TEnergy >::
PerformSwap( LabelType aLabel, LabelType bLabel ) {
	EnergyValType afterSwapEnergy = 0.0;
	this->m_ActiveSites.clear();
    LabelType* lBuffer = this->m_MRFEnergy->GetOutputPointer();
	SiteVector idHolder = this->m_MRFEnergy->GetIdHolder();
	LabelType cLabel;

	// Build list of active sites based on alpha and current labeling
	typename SiteVector::const_iterator it;
	typename SiteVector::const_iterator end = idHolder.end();

	//for (it = idHolder.begin(); it!=end; it++ ) {
	for ( size_t i = 0, lookupSite = 0; i < idHolder.size(); i++) {
		size_t offset = idHolder[i];
		cLabel = *( lBuffer + offset );

		assert( cLabel!= 0 );

		if ( cLabel == aLabel || cLabel == bLabel ) {
			this->m_ActiveSites.push_back( offset );
			this->m_LookupSiteVar[ offset ] = lookupSite++;
        }
    }

	size_t size = this->m_ActiveSites.size();

	this->m_BeforeMoveEnergy = 0.0;

	if ( size > 0 ) {
		GraphType *e = new GraphType(size, 6 ); // TODO fix number of neighbors
		e->add_variable(size);

        SetUpDataCostsSwap( aLabel, bLabel, e );
        SetUpSmoothCostsSwap( aLabel, bLabel, e);

        afterSwapEnergy = e->minimize();

        if( afterSwapEnergy < this->m_BeforeMoveEnergy ) {    	// Relabeling
        	size_t changed = 0;
        	end = this->m_ActiveSites.end();

	        for ( it = this->m_ActiveSites.begin(); it!=end ; ++it ) {
	        	size_t offset = *it;
	        	*( lBuffer + offset ) = ( e->get_var( this->m_LookupSiteVar[offset] )==0 ) ?
	        		                    aLabel : bLabel;
	        	changed++;
	        }
	        if  (changed>0) this->m_MRFEnergy->LabelingModified();
        }

        this->ResetLookupSiteVar();
        delete e;
    }
};

template < class TEnergy >
void
SwapOptimizer< TEnergy >::
SetUpDataCostsSwap( LabelType aLabel, LabelType bLabel, GraphType *e ) {
	EnergyValType alphaT, betaT;

	typename SiteVector::const_iterator end = (this->m_ActiveSites).end();
	size_t i = 0;
	for (typename SiteVector::const_iterator it = (this->m_ActiveSites).begin(); it!=end; it++, i++) {
    	alphaT = this->m_MRFEnergy->GetSampleEnergyForLabel(aLabel, *it);
    	betaT  = this->m_MRFEnergy->GetSampleEnergyForLabel(bLabel, *it);

    	this->AddtLink(e, i, alphaT,betaT);
    }
}

template < class TEnergy >
void
SwapOptimizer< TEnergy >::
SetUpSmoothCostsSwap( LabelType aLabel, LabelType bLabel, GraphType *e ) {
	// Prepare some needed variables
	size_t offset;
	size_t nOffset = -1;

	size_t size = this->m_ActiveSites.size();
	LabelType* lBuffer = this->m_Output->GetBufferPointer();

	typename LabelImageType::SizeType maxIndex = this->m_Output->GetLargestPossibleRegion().GetSize();
	typename LabelImageType::IndexType siteIndex;
	typename LabelImageType::IndexType nSiteIndex;

	LabelType cLabel;
	LabelType nLabel;

	for(int i = size-1; i>=0; --i ) {
	    // Get pixel offset index from the lookup-variable array
		offset = this->m_ActiveSites[i];

		// Compute current's pixel index and current label
	    siteIndex = this->m_Output->ComputeIndex( offset );
		cLabel = *( lBuffer + offset );

		// Search in direct neighborhood and create the energy graph
	    for (size_t dim = 0; dim < ImageDimension; dim++) {
	    	for( size_t inc = 0; inc < 2; inc++ ) {
				nSiteIndex = siteIndex;
				nSiteIndex[dim] += inc*2 - 1;

				if ( nSiteIndex[dim] >= 0 && nSiteIndex[dim] < maxIndex[dim] ) {
					nOffset = this->m_Output->ComputeOffset(nSiteIndex);
					nLabel = *( lBuffer + nOffset );

					int siteVarId = this->m_LookupSiteVar[nOffset];

					bool isInLookup = (siteVarId != -1 );

					if (!isInLookup && aLabel!=nLabel && bLabel!=nLabel) {
						EnergyValType E01 = this->m_MRFEnergy->GetNeighborCost(aLabel, nLabel );
						EnergyValType E10 = this->m_MRFEnergy->GetNeighborCost(bLabel, nLabel );
						this->AddtLink(e,i, E01, E10 );
					} else if (isInLookup && (nOffset < offset) ) {
						EnergyValType E00 = this->m_MRFEnergy->GetNeighborCost(aLabel, aLabel);
						EnergyValType E01 = this->m_MRFEnergy->GetNeighborCost(aLabel, bLabel);
						EnergyValType E10 = this->m_MRFEnergy->GetNeighborCost(bLabel, aLabel);
						EnergyValType E11 = this->m_MRFEnergy->GetNeighborCost(bLabel, bLabel);

						this->AddnLink(e,i,siteVarId,E00, E01, E10, E11 );
					}
				}
	    	}
	    }
	}
};

#endif /* MRFSWAPOPTIMIZER_TXX_ */
