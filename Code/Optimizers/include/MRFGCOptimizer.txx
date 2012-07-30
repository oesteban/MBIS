// --------------------------------------------------------------------------------------
// File:    MRFGCOptimizer.txx
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

#ifndef MRFGCOPTIMIZER_TXX_
#define MRFGCOPTIMIZER_TXX_

#include "MRFGCOptimizer.h"

#include "itkNumericTraits.h"


template < class TEnergy >
MRFGCOptimizer<TEnergy>::
MRFGCOptimizer(): m_BeforeMoveEnergy(0.0), m_MaxAlgorithmIterations(5), m_NumberOfClasses(0)
{

}


template < class TEnergy >
void MRFGCOptimizer<TEnergy>::
StartOptimization() {
	if ( this->m_MRFEnergy.IsNull() ) { // TODO throw exception
		std::cerr << "StartOptimization: EnergyFunction not provided" << std::endl;
		exit(1);
	}

	this->m_Output = this->m_MRFEnergy->GetOutput();
	this->m_LookupSiteVar.resize( this->m_Output->GetLargestPossibleRegion().GetNumberOfPixels() );
	this->ResetLookupSiteVar();

	this->m_NumberOfClasses = this->m_MRFEnergy->GetNumberOfClasses();
	this->m_LabelTable.SetSize( this->m_NumberOfClasses );
	for (size_t i = 1; i <= this->m_NumberOfClasses; i++)
		this->m_LabelTable[i-1] = i;

	EnergyValType Ed, Es, Etotal, Ed_old, Es_old, Etotal_old;
	Ed_old = Ed = this->m_MRFEnergy->GetTotalDataEnergy();
	Es_old = Es = this->m_MRFEnergy->GetTotalSmoothnessEnergy();
	Etotal_old = Etotal = Ed + Es;
    int noChange = 0;

	std::cout << "\t\tInit: ";
	std::cout << "Et = " << std::setprecision(15) << Etotal << " (Ed=" << Ed << ", Es=" << Es << ")." << std::endl << std::setprecision(4);

	itk::TimeProbe clock;

    for (size_t iter = 0; iter < m_MaxAlgorithmIterations; iter++) {
    	clock.Start();
    	optimize();
    	clock.Stop();
		Ed = this->m_MRFEnergy->GetTotalDataEnergy();
		Es = this->m_MRFEnergy->GetTotalSmoothnessEnergy();
		Etotal = Ed + Es;

		std::cout << "\t\t[" << std::setfill('0') << std::setw(2) << (iter+1) << "]: " << std::setfill(' ');
		std::cout << "Et = " << std::setprecision(15) << Etotal << " (Ed=" << Ed << ", Es=" << Es << ")." << std::setprecision(4);
		if (Etotal > Etotal_old)  std::cout << " Warning: energy is increasing!.";
		std::cout << std::endl;

        if ((Etotal == Etotal_old) && (Ed == Ed_old) && (Es==Es_old)) {
        	noChange++;
        	if ( noChange >= 2 ) break;
        } else
        	noChange = 0;

		Etotal_old = Etotal;
        Ed_old = Ed;
        Es_old = Es;
    }

    std::cout << "\t\t      Segmentation finished. Time: mean=" << clock.GetMean() << ", total=" << clock.GetTotal()<< " (s)" << std::endl;

}

template < class TEnergy >
void MRFGCOptimizer<TEnergy>::
optimize() {
	// Boykov et al. (2001) implementation of expansion (Fig.3 paper)
	EnergyValType Ed = this->m_MRFEnergy->GetTotalDataEnergy();
	EnergyValType Es = this->m_MRFEnergy->GetTotalSmoothnessEnergy();
	EnergyValType new_energy = Ed + Es;
	EnergyValType old_energy = itk::NumericTraits<EnergyValType>::max();

	size_t current_it = 1;

	while ( old_energy > new_energy && current_it <= this->m_MaxAlgorithmIterations ) {
		old_energy = new_energy;
		this->ScrambleLabelTable();

		// 3. For each label alpha, do alpha-expansion
		this->PerformMove();

		Ed = this->m_MRFEnergy->GetTotalDataEnergy();
		Es = this->m_MRFEnergy->GetTotalSmoothnessEnergy();
		new_energy = Ed + Es;
		current_it++;
	}
}




template < class TEnergy >
inline void
MRFGCOptimizer<TEnergy>::
AddtLink( GraphType *e, size_t i, EnergyValType e0, EnergyValType e1) {
	if ( e0 > GCO_MAX_ENERGYTERM || e1 > GCO_MAX_ENERGYTERM ) {
		std::cerr << "Cost term was larger than GCO_MAX_ENERGYTERM; danger of integer overflow." << std::endl;
		exit(1);
	}
	if ( e0 < 0 || e1 < 0 ) {
		std::cerr << "Cost term was negative" << std::endl;
	}
	m_BeforeMoveEnergy+=e1;

	e->add_term1( i, e0, e1 );
}

template < class TEnergy >
inline void
MRFGCOptimizer<TEnergy>::
AddnLink( GraphType *e, size_t i, size_t j, EnergyValType e00, EnergyValType e01, EnergyValType e10, EnergyValType e11 )
{
	if ( e00 > GCO_MAX_ENERGYTERM || e11 > GCO_MAX_ENERGYTERM || e01 > GCO_MAX_ENERGYTERM || e10 > GCO_MAX_ENERGYTERM ) {
		std::cerr << "Smooth cost term was larger than GCO_MAX_ENERGYTERM; danger of integer overflow." << std::endl;
		exit(1);
	}
	if ( e00 < 0 || e11 < 0 || e01 < 0 || e10 < 0 ) {
		std::cerr << "Smooth cost term was negative" << std::endl;
	}
	if ( e00+e11 > e01+e10 ) {
		std::cerr<< "Non-submodular expansion term detected; smooth costs must be a metric for expansion" << std::endl;
		exit(1);
	}

	m_BeforeMoveEnergy += e11;

	e->add_term2( i,j,e00,e01,e10,e11 );
};

template < class TEnergy >
void MRFGCOptimizer<TEnergy>::
ScrambleLabelTable() {
   size_t num_times = this->m_NumberOfClasses*2;
   LabelType r1,r2,temp;

   for (size_t cnt = 0; cnt < num_times; cnt++ ) {
      r1 = rand()%this->m_NumberOfClasses;
      r2 = rand()%this->m_NumberOfClasses;

      temp                   = this->m_LabelTable[r1];
      this->m_LabelTable[r1] = this->m_LabelTable[r2];
      this->m_LabelTable[r2] = temp;
   }
}


#endif /* MRFGCOPTIMIZER_TXX_ */
