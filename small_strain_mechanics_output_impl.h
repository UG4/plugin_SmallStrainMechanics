/*
 * small_strain_mechanics_output_impl.h
 *
 *  Created on: 22.03.2013
 *      Author: raphaelprohl
 */

#ifndef SMALL_STRAIN_MECHANICS_OUTPUT_IMPL_H_
#define SMALL_STRAIN_MECHANICS_OUTPUT_IMPL_H_

// other ug4 modules
#include "common/common.h"

#include "lib_disc/common/local_algebra.h"
#include "lib_grid/geometric_objects/geometric_objects.h"

// module intern headers
#include "small_strain_mechanics_output.h"

namespace ug {
namespace SmallStrainMechanics{

template <typename TDomain, typename TGridFunction>
SmallStrainMechOutput<TDomain,TGridFunction>::SmallStrainMechOutput(
		SmartPtr<SmallStrainMechanicsElemDisc<TDomain> > spElemDisc)
: m_spElemDisc(spElemDisc)
{}

template <typename TDomain, typename TGridFunction>
void
SmallStrainMechOutput<TDomain,TGridFunction>::normal_stresses_strains(TGridFunction& devSigma,
		TGridFunction& sigma, TGridFunction& epsilon, const TGridFunction& u)
{
	static const int dim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object geometric_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;

	// 	local indices and local algebra
	LocalIndices indU, indEps, indSig, indDevSig;
	LocalVector locU, locSig, locEps, locDevSig;

	const_iterator iter = u.template begin<geometric_base_object>();
	const_iterator end = u.template end<geometric_base_object>();

	//	loop over all elements
	for(;iter != end; ++iter)
	{
		//	get element
		geometric_base_object* elem = *iter;

		// 	get global indices
		u.indices(elem, indU); epsilon.indices(elem, indEps);
		sigma.indices(elem, indSig); devSigma.indices(elem, indDevSig);

		// 	adapt local algebra
		locU.resize(indU); locEps.resize(indEps);
		locSig.resize(indSig); locDevSig.resize(indDevSig);

		//	reset contribution of this element
		locDevSig = 0.0; locSig = 0.0; locEps = 0.0;

		//	local vector extract -> locU
		GetLocalVector(locU, u);

		//	call method of base class SmallStrainMechanicsElemDisc
		m_spElemDisc->normal_stress_strain_loc(locDevSig, locSig, locEps, elem, locU);

		// 	send local to global
		//	sigma (cauchy-stress tensor), epsilon (linearized strain tensor), devSigma (deviatoric part of Sigma)
		AddLocalVector(devSigma, locDevSig);
		AddLocalVector(sigma, locSig);
		AddLocalVector(epsilon, locEps);
	}

}

} // namespace SmallStrainMechanics
} // namespace ug

#endif /* SMALL_STRAIN_MECHANICS_OUTPUT_IMPL_H_ */
