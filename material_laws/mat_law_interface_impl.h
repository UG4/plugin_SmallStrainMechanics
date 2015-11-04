/*
 * mat_law_interface_impl.h
 *
 *  Created on: 03.02.2014
 *      Author: raphaelprohl
 */

#ifndef MAT_LAW_INTERFACE_IMPL_H_
#define MAT_LAW_INTERFACE_IMPL_H_

#include "mat_law_interface.h"

namespace ug{
namespace SmallStrainMechanics{

template<typename TDomain>
template<typename TFEGeom>
void
IMaterialLaw<TDomain>::
DisplacementGradient(MathMatrix<dim, dim>& GradU, const size_t ip,
		const TFEGeom& geo, const LocalVector& u)
{
	//	loop shape-functions at one integration point ip in order
	//	to compute local_grad(ip,i): \frac{\partial N_i}{\eps_ip}
	//	and global_grad(ip,i): \frac{\partial N_i}{\X_ip}

	for (size_t i = 0; i < (size_t) dim; ++i)
		for (size_t J = 0; J < (size_t) dim; ++J)
		{
			GradU[i][J] = 0.0;

			//	compute GradU: displacementGradient
			for (size_t a = 0; a < geo.num_sh(); ++a)
				GradU[i][J] += geo.global_grad(ip, a)[J] * u(i, a);
		}
}

}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#endif /* MAT_LAW_INTERFACE_IMPL_H_ */
