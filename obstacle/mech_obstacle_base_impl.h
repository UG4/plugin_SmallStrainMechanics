/*
 * Copyright (c) 2026:  KAUST
 * Authors: R. Krause, D. Logashenko
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

namespace ug{
namespace SmallStrainMechanics{

/* Class BoxVertexConstraintManager */

/**
 * Marks the Signorini boundary DoFs and computes the normal vectors (for one element type)
 */
template <typename TDomain, typename TAlgebra>
template <typename TBaseElem>
void SignoriniConstraint<TDomain, TAlgebra>::mark_signorini_bnd_elem
(
	const DoFDistribution& dd, ///< DoF distribution of the solution
	FunctionGroup& fctGrp, ///< function group of the displacement
	matrix_type& transform, ///< the transformation (rotation) matrix
	vector_type& bnd_gf, ///< distribution of the bound (and the flags)
	number time ///< the time argument
)
{
	typedef typename DoFDistribution::traits<TBaseElem>::const_iterator e_iter_t;
	
	if(fctGrp.size() != dim)
		UG_THROW("SignoriniConstraint: Number of the components is not equal to the dimensionality");
	
	const LFEID& lfeID = dd.local_finite_element_id(fctGrp[0]);
	for(size_t i = 1; i < fctGrp.size(); i++)
		if(dd.local_finite_element_id(fctGrp[i]) != lfeID)
			UG_THROW("SignoriniConstraint: Components of the grid function should be of the same type!");
	
	std::vector<DoFIndex> multInd[dim];
	std::vector<position_type> vPos;
	std::vector<number> bound(1);
	std::vector<MathVector<dim> > nrm_dir(1);

	for(const auto& c_pair : constr_maps())
	{
		int si = c_pair.first;
		SmartPtr<SignoriniObstacle<domain_type> > constr = c_pair.second;
		
		e_iter_t e_end = dd.template end<TBaseElem>(si);
		for(e_iter_t e_iter = dd.template begin<TBaseElem>(si); e_iter != e_end; ++e_iter)
		{
			TBaseElem* elem = *e_iter;
			
		//	get dof positions
			InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);
			// We assume that the DoF positions (dp's) in vPos are stored in the same order
			// as the indices returned by dd->inner_dof_indices for every function
			// in fctGrp. The same assumption is used in lagrange_dirichlet_boundary_impl.h.
			
			size_t n_dp = vPos.size(); // number of the associated DoFs
			
		//	get the multiindices for every function
			for(size_t f = 0; f < dim; f++)
			{
				dd.inner_dof_indices(elem, fctGrp[f], multInd[f]);
				if(multInd[f].size() != n_dp)
					UG_THROW("SignoriniConstraint: Index number mismatch.");
			}

		//	get the bounds and the directions
			bound.resize(vPos.size());
			nrm_dir.resize(vPos.size());
			(* (constr->value())) (&(bound[0]), &(vPos[0]), time, si, n_dp);
			(* (constr->normal())) (&(nrm_dir[0]), &(vPos[0]), time, si, n_dp);
		
		//	set the flags, save the directions and the bounds
			//	We mark only the first component in every DoF position:
			//	After the rotation, only the first component is restricted.
			//	For the same reason, we also store the bound only for the first component.
			for(size_t dp = 0; dp < n_dp; dp++)
			{
				DoFFlagManager<algebra_type>::set_flag
					(bnd_gf, multInd[0][dp][0],
						m_bnd_flag[0], m_bnd_flag[1], multInd[0][dp][1]);
				
				DoFRef(bnd_gf, multInd[0][dp]) = bound[dp];
			}
			
		//	compute the transformation
			MathMatrix<dim,dim> dofTransform;
			
			for(size_t dp = 0; dp < n_dp; dp++)
			{
				MathVector<dim> dir = nrm_dir[dp];
				dir[0] += 1;
				MatHouseholder(dofTransform, dir);
				
				for(size_t f_i = 0; f_i < dim; f_i++)
					for(size_t f_j = 0; f_j < dim; f_j++)
						DoFRef(transform, multInd[f_i][dp], multInd[f_j][dp])
							= dofTransform(f_i, f_j);
			}
		}
	}
}

/**
 * Marks the Signorini boundary DoFs and computes the normal vectors
 */
template <typename TDomain, typename TAlgebra>
void SignoriniConstraint<TDomain, TAlgebra>::mark_signorini_bnd
(
	const DoFDistribution& dd, ///< DoF distribution of the solution
	FunctionGroup& fctGrp, ///< function group of the displacement
	matrix_type& transform, ///< the transformation (rotation) matrix
	vector_type& bnd_gf, ///< distribution of the bound (and the flags)
	number time ///< the time argument
)
{
	if(fctGrp.size() != dim)
		UG_THROW("SignoriniConstraint: Number of the functions is not equal to the dimensionality.");
	
//	check if the flags are present; if not, create them
	DoFFlagManager<algebra_type>::get()->create_flags(bnd_gf);
	for(size_t i = 0; i < bnd_gf.size(); i++)
		DoFFlagManager<algebra_type>::set_flag_block
			(bnd_gf, i, m_bnd_flag[0], m_bnd_flag[1], 0);
	
//	initialize the transform and the bound
	transform.set(1.0); // set to the identity matrix
	bnd_gf.set(0.0); // set to zero to initialize everywhere (it is used only where the flags are set)
	
//	set the transform and the bound at the selected boundary
	if(dd.max_dofs(VERTEX))
		mark_signorini_bnd_elem<RegularVertex>(dd, fctGrp, transform, bnd_gf, time);
	if(dd.max_dofs(EDGE))
		mark_signorini_bnd_elem<Edge>(dd, fctGrp, transform, bnd_gf, time);
	if(dd.max_dofs(FACE))
		mark_signorini_bnd_elem<Face>(dd, fctGrp, transform, bnd_gf, time);
	if(dd.max_dofs(VOLUME))
		mark_signorini_bnd_elem<Volume>(dd, fctGrp, transform, bnd_gf, time);
}

} //end of namespace SmallStrainMechanics
} //end of namespace ug

/* End of File */