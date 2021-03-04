/*
 * Copyright (c) 2014-2020:  G-CSC, Goethe University Frankfurt
 * Author: Lukas Larisch
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

#ifndef ANISO_REF_
#define ANISO_REF_

#include "lib_grid/algorithms/refinement_mark_util.h"
#ifdef UG_PARALLEL
#include "lib_grid/parallelization/util/attachment_operations.hpp"
#include "lib_grid/parallelization/parallel_refinement/parallel_refinement.h"
#endif

#include "common/math/misc/math_util.h"

namespace ug{
namespace SmallStrainMechanics{

#ifdef UG_DIM_3

template <class TDomain>
number MaxEdgeLength(TDomain& dom, MathVector<TDomain::dim>& ndir)
{
	UG_ASSERT(dom.grid().get() == refiner.get_associated_grid(),
			  "Grids in domain and in refiner have to match!");

	UG_ASSERT(TDomain::dim == 3, "dimension must be 3!");

	typedef typename domain_traits<TDomain::dim>::element_type elem_t;
	typedef typename MultiGrid::traits<elem_t>::iterator iter_t;

	typename TDomain::position_accessor_type aaPos = dom.position_accessor();
	MultiGrid& mg = *dom.grid();
	MultiGrid::edge_traits::secure_container edges;

	typedef MathVector<TDomain::dim> vector_t;

	number maxSq = 0;

	for(iter_t e_iter = mg.begin<elem_t>(); e_iter != mg.end<elem_t>(); ++e_iter){
		elem_t* elem = *e_iter;
		if(mg.has_children(elem))
			continue;

		mg.associated_elements(edges, elem);

		for(size_t i = 0; i < edges.size(); ++i){
			vector_t edgeDir;
			VecSubtract(edgeDir, aaPos[edges[i]->vertex(0)], aaPos[edges[i]->vertex(1)]);
			VecNormalize(edgeDir, edgeDir);
			number dot = VecDot(edgeDir, ndir);

			if((dot + SMALL >= 1) || (dot - SMALL <= -1)){
				maxSq = max(maxSq, EdgeLengthSq(edges[i], aaPos));
			}
		}
	}

	return sqrt(maxSq);
}



template <class TDomain>
void MarkAnisotropic_Along_Normal(TDomain& dom, IRefiner& refiner, MathVector<TDomain::dim>& ndir)
{
	UG_ASSERT(dom.grid().get() == refiner.get_associated_grid(),
			  "Grids in domain and in refiner have to match!");

	UG_ASSERT(TDomain::dim == 3, "dimension must be 3!");

	typedef typename domain_traits<TDomain::dim>::element_type elem_t;
	typedef typename MultiGrid::traits<elem_t>::iterator iter_t;

	typename TDomain::position_accessor_type aaPos = dom.position_accessor();
	MultiGrid& mg = *dom.grid();
	MultiGrid::edge_traits::secure_container edges;

	typedef MathVector<TDomain::dim> vector_t;

	for(iter_t e_iter = mg.begin<elem_t>(); e_iter != mg.end<elem_t>(); ++e_iter){
		elem_t* elem = *e_iter;
		if(mg.has_children(elem))
			continue;

	//	we'll mark all volumes as anisotropic, since we currently have to use
	//	copy elements during anisotropic refinement.
		refiner.mark(elem, RM_ANISOTROPIC);

		mg.associated_elements(edges, elem);

		for(size_t i = 0; i < edges.size(); ++i){
			vector_t edgeDir;
			VecSubtract(edgeDir, aaPos[edges[i]->vertex(0)], aaPos[edges[i]->vertex(1)]);
			VecNormalize(edgeDir, edgeDir);
			number dot = VecDot(edgeDir, ndir);

			if((dot + SMALL >= 1) || (dot - SMALL <= -1)){
				refiner.mark(edges[i]);
			}
		}
	}
}


template <class TDomain>
void MarkAnisotropic_Longest_Scaled_Normal(TDomain& dom, IRefiner& refiner, number Ex, number Ey, number Ez)
{
	UG_ASSERT(dom.grid().get() == refiner.get_associated_grid(),
			  "Grids in domain and in refiner have to match!");

	UG_ASSERT(TDomain::dim == 3, "dimension must be 3!");

	typedef MathVector<TDomain::dim> vector_t;

	vector_t ndirX; //x tang
	vector_t ndirY; //y long
	vector_t ndirZ; //z radi

	ndirX[0] = 1; ndirX[1] = 0; ndirX[2] = 0;
	ndirY[0] = 0; ndirY[1] = 1; ndirY[2] = 0;
	ndirZ[0] = 0; ndirZ[1] = 0; ndirZ[2] = 1;

	//length of edges in directions x, y, z
	number maxX = MaxEdgeLength(dom, ndirX);
	number maxY = MaxEdgeLength(dom, ndirY);
	number maxZ = MaxEdgeLength(dom, ndirZ);

	number scaleXsqrt = sqrt(maxX/Ex);
	number scaleYsqrt = sqrt(maxY/Ey);
	number scaleZsqrt = sqrt(maxZ/Ez);

	number scaleXYZmin = min(scaleXsqrt, min(scaleYsqrt, scaleZsqrt));

	scaleXsqrt /= scaleXYZmin;
	scaleYsqrt /= scaleXYZmin;
	scaleZsqrt /= scaleXYZmin;

	std::cout << "scale normed x: " << scaleXsqrt << std::endl;
	std::cout << "scale normed y: " << scaleYsqrt << std::endl;
	std::cout << "scale normed z: " << scaleZsqrt << std::endl;

	if(scaleXsqrt > scaleYsqrt){
		if(scaleXsqrt > scaleZsqrt){
			UG_LOG("[MarkAnisotropic_Longest_Scaled_Normal] refine in X direction");
			MarkAnisotropic_Along_Normal(dom, refiner, ndirX);
		}
		else{
			UG_LOG("[MarkAnisotropic_Longest_Scaled_Normal] refine in Z direction");
			MarkAnisotropic_Along_Normal(dom, refiner, ndirZ);
		}

	}
	else{
		if(scaleZsqrt > scaleYsqrt){
			UG_LOG("[MarkAnisotropic_Longest_Scaled_Normal] refine in Z direction");
			MarkAnisotropic_Along_Normal(dom, refiner, ndirZ);
		}
		else{
			UG_LOG("[MarkAnisotropic_Longest_Scaled_Normal] refine in Y direction");
			MarkAnisotropic_Along_Normal(dom, refiner, ndirY);
		}
	}
}

#endif


}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#endif /* ANISO_REF_ */
