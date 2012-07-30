// --------------------------------------------------------------------------------------
// File:          MRFExpansionOptimizer.txx
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

#ifndef MRFEXPANSIONOPTIMIZER_TXX_
#define MRFEXPANSIONOPTIMIZER_TXX_

#include "MRFExpansionOptimizer.h"

template < class TEnergy >
ExpansionOptimizer< TEnergy >::
ExpansionOptimizer() {

};

template < class TEnergy >
void
ExpansionOptimizer< TEnergy >::
PerformAlphaExpansion( LabelType aLabel ) {
	EnergyValType afterExpansionEnergy = 0.0;
	this->m_ActiveSites.clear();

	LabelType * lBuffer = this->m_MRFEnergy->GetOutputPointer();
	size_t offset;

	SiteVector idHolder = this->m_MRFEnergy->GetIdHolder();

	// Build list of active sites based on alpha and current labeling
	typename SiteVector::const_iterator it;
	typename SiteVector::const_iterator end = idHolder.end();
	LabelType cLabel;

	//for (it = idHolder.begin(); it!=end; it++ ) {
	for ( size_t i = 0, lookupSite = 0; i < idHolder.size(); i++) {
		size_t offset = idHolder[i];
		cLabel = *( lBuffer + offset );
		assert( cLabel != 0 );

		if ( cLabel != aLabel ) {
			this->m_ActiveSites.push_back( offset );
			this->m_LookupSiteVar[offset] = lookupSite;
			lookupSite++;
        }
    }

	size_t size = this->m_ActiveSites.size();

	// Initialize reverse-lookup so that non-active neighbors can be identified
	// while constructing the graph
	//for ( size_t i=0; i < size; i++)
	//	this->m_LookupSiteVar[this->m_ActiveSites[i]]=i;

	this->m_BeforeMoveEnergy = 0.0;

	if ( size > 0 ) {
		GraphType *e = new GraphType(size+this->m_NumberOfClasses,size + this->m_NumberOfClasses + 6 ); // TODO fix number of neighbors
		e->add_variable(size);

        SetUpDataCostsExpansion( aLabel,e );
        SetUpSmoothCostsExpansion( aLabel,e );

        EnergyValType alphaCorrection = 0.0; // TODO setupLabelCostsExpansion method

        afterExpansionEnergy = e->minimize() + alphaCorrection;
        typename LabelImageType::IndexType siteIdx;
        if( afterExpansionEnergy < this->m_BeforeMoveEnergy ) {
        	size_t changed = 0;

        	// Relabeling
        	end = this->m_ActiveSites.end();
	        for ( it = this->m_ActiveSites.begin(), size=0; it!=end ; it++ ) {
	        	if ( e->get_var( this->m_LookupSiteVar[*it] )==0 ) {
	        		*( lBuffer + *it ) = aLabel;
	        		changed++;
		    	}
	        }
	        if  (changed>0) this->m_MRFEnergy->LabelingModified();
        }
        this->ResetLookupSiteVar();

        delete e;
    }
}

template < class TEnergy >
void
ExpansionOptimizer< TEnergy >::
SetUpDataCostsExpansion( LabelType aLabel, GraphType *e ) {
	LabelType cLabel;
	EnergyValType alphaT, labelT;

	typename SiteVector::const_iterator end = (this->m_ActiveSites).end();
	size_t i = 0;
	for (typename SiteVector::const_iterator it = (this->m_ActiveSites).begin(); it!=end; it++, i++) {
    	alphaT = this->m_MRFEnergy->GetSampleEnergyForLabel(aLabel, *it);
    	labelT = this->m_MRFEnergy->GetSampleEnergy(*it);
    	this->AddtLink(e, i, alphaT,labelT);
    }
}

template < class TEnergy >
void
ExpansionOptimizer< TEnergy >::
SetUpSmoothCostsExpansion(LabelType aLabel,GraphType *e ) {
	// Prepare some needed variables
    size_t offset;
	size_t nOffset = -1;

	size_t size = this->m_ActiveSites.size();
    LabelType* lBuffer = this->m_MRFEnergy->GetOutputPointer();

    typename LabelImageType::SizeType maxIndex = this->m_Output->GetLargestPossibleRegion().GetSize();
	typename LabelImageType::IndexType siteIndex;
	typename LabelImageType::IndexType nSiteIndex;

	LabelType cLabel;
	LabelType nLabel;

    for(int i = size-1; i>=0; --i ) {
        // Get pixel offset index from the lookup-variable array
		offset = this->m_ActiveSites[i];

		// Compute current's pixel index
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

					if (!isInLookup && nLabel == aLabel ) {
						EnergyValType E01 = this->m_MRFEnergy->GetNeighborCost(aLabel, nLabel );
						EnergyValType E10 = this->m_MRFEnergy->GetNeighborCost(cLabel, nLabel );
						this->AddtLink(e,i, E01, E10 );
					} else if (isInLookup && (nOffset < offset) ) {
						EnergyValType E00 = this->m_MRFEnergy->GetNeighborCost(aLabel, aLabel );
						EnergyValType E01 = this->m_MRFEnergy->GetNeighborCost(aLabel, nLabel );
						EnergyValType E10 = this->m_MRFEnergy->GetNeighborCost(cLabel, aLabel );
						EnergyValType E11 = this->m_MRFEnergy->GetNeighborCost(cLabel, nLabel );

						this->AddnLink(e,i,siteVarId,E00, E01, E10, E11 );
					}
				}
        	}
        }
    }
}


#endif /* MRFEXPANSIONOPTIMIZER_TXX_ */
