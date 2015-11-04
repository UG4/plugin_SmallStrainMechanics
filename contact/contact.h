/*
 * contact.h
 *
 *  Created on: 06.05.2013
 *      Author: raphaelprohl
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
