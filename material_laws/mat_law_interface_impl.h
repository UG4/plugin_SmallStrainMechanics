/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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
