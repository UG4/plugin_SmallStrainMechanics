/*
 * adaptive_util.h
 *
 *  Created on: 23.03.2015
 *      Author: raphaelprohl
 */

#ifndef ADAPTIVE_UTIL_H_
#define ADAPTIVE_UTIL_H_

#include "lib_grid/algorithms/refinement/refiner_interface.h"

namespace ug{
namespace SmallStrainMechanics{

template<typename TDomain, typename TElem>
void plastic_ip_elem(bool& bPlasticIPs, TElem* elem,
		const LocalVector& locU, SmartPtr<TDomain> dom,
		SmallStrainMechanicsElemDisc<TDomain>& elemDisc);

template <typename TDomain, typename TAlgebra>
void MarkForAdaption_PlasticElem(IRefiner& refiner,
     GridFunction<TDomain, TAlgebra>& u,
     SmallStrainMechanicsElemDisc<TDomain>& elemDisc);

} //end of namespace SmallStrainMechanics
} //end of namespace ug

#include "adaptive_util_impl.h"

#endif /* ADAPTIVE_UTIL_H_ */
