/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef CONTACT_H_
#define CONTACT_H_

#include "bridge/bridge.h"

// other ug4 modules
#include "common/common.h"
#include "lib_algebra/active_set/lagrange_multiplier_disc_interface.h"

// intern headers
#include "../small_strain_mech.h"

using namespace std;

namespace ug{
namespace SmallStrainMechanics{

/// Contact conditions for small strain mechanical applications
/**
 * Begin of the considering contact conditions in a structural mechanical model
 * for small strains. This implementation is based on active sets.
 *
 * Still work in progress. A first introductive work for active sets and Lagrange
 * multiplier can be found in:
 *
 * References:
 * <ul>
 * <li> Corinna Hager und Babara I. Wohlmuth: Hindernis- und Kontaktprobleme
 * <ul>
 */

template <typename TDomain, typename TGridFunction>
class ContactSmallStrainMechanics
	:  public ILagrangeMultiplierDisc<TDomain, TGridFunction>
{
	private:
	///	Base class type
		typedef ILagrangeMultiplierDisc<TDomain, TGridFunction> base_type;

	///	own type
		typedef ContactSmallStrainMechanics<TDomain, TGridFunction> this_type;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	public:
		ContactSmallStrainMechanics(
				SmartPtr<SmallStrainMechanicsElemDisc<TDomain> > spLinElastPlast);

	/// Virtual destructor
		virtual ~ContactSmallStrainMechanics() {}

		virtual void lagrange_multiplier(TGridFunction& lagMult, const TGridFunction& u,
				vector<DoFIndex> vActiveSet, vector<int> vActiveSubsets);

		template <typename TSideElem, typename TIterator>
		void contact_forces_elem(TIterator iterBegin,
				TIterator iterEnd, const TGridFunction& u, TGridFunction& contactForce,
				vector<DoFIndex> vActiveSetGlob);

	private:
	///	pointer to the element discretization
		SmartPtr<SmallStrainMechanicsElemDisc<TDomain> >  m_spMechElemDisc;
};

}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#include "contact_impl.h"

#endif /* CONTACT_H_ */
