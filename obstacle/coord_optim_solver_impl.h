/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko, Rolf Krause
 * Based on: linear_solver.h by Andreas Vogel and lu.h by Martin Rupp
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

/**
 * Compute the correction \f$ c \gets L^{-1} d \f$ by the constrained Gauss-Seidel
 * method
 */
template <typename TAlgebra>
bool ObstacleGaussSeidel<TAlgebra>::compute_correction
(
	vector_type& c,	///< the correction to compute
	const vector_type& d,	///< the provided defect
	const vector_type& u	///< the solution (previous iterate)
)
{
	typedef typename matrix_type::value_type matrix_block;
	typedef typename vector_type::value_type vector_block;
	typedef typename matrix_type::const_row_iterator const_row_it;
	
	bool constrain = false;
	
	if(m_spLowerBnd.valid() || m_spUpperBnd.valid())
	{
		constrain = true;
		
		if(m_obst_bnd_flag_name.empty())
			UG_THROW("Bounds specified but no flag name is given!");
		if(m_spLowerBnd.valid() && m_spLowerBnd->size() != u.size())
			UG_THROW("Size of the lower bound vector is not equal to the size of the solution vector!");
		if(m_spUpperBnd.valid() && m_spUpperBnd->size() != u.size())
			UG_THROW("Size of the upper bound vector is not equal to the size of the solution vector!");
		if(! u.flag_set_present(m_obst_bnd_flag[0]))
			UG_THROW("Flag '" << m_obst_bnd_flag_name << "' is not present in the solution");
	}
	
	const size_t vec_len = c.size();
	
	vector_block s, t;

	for(size_t i = 0; i < vec_len; i++)
	{
		s = d[i];

		const const_row_it rowEnd = m_pA->end_row(i);
		const_row_it it = m_pA->begin_row(i);
		for(; it != rowEnd && it.index() < i; ++it) // for j < i (where j := it)
			// s -= A[i,j] * c[j];
			MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()]);
		
		if(it.index() != i)
			return false; // no diagonal block found!!!
		matrix_block A = it.value();

		// t = s / A[i,i]
		InverseMatMult(t, 1.0, A, s);
		
		// check if we need to constrain the correction
		flag_unit_type flags;
		if(constrain)
			flags = flagBlockFrame & flag_mng_type::flag_block
				(u, i, m_obst_bnd_flag[0], m_obst_bnd_flag[1]);
		else
			flags = 0;
		
		if(flags)
		{
			if(m_spLowerBnd.valid()) // constrain the correction from below
			{
				vector_block diff = (* m_spLowerBnd)[i] - u[i];
				for(int cmp = 0; cmp < blockSize; cmp++)
				{
					if(((flags >> cmp) & 1)
						&& BlockRef(t, cmp) > BlockRef(diff, cmp))
					{
						// set the "Dirichlet row" and recompute the correction
						for(size_t k = 0; k < blockSize; k++) BlockRef(A, i, k) = 0;
						BlockRef(A, i, i) = 1;
						BlockRef(s, i) = BlockRef(diff, i);
						InverseMatMult(t, 1.0, A, s);
					}
				}
			}
			
			if(m_spUpperBnd.valid()) // constrain the correction from above
			{
				vector_block diff = (* m_spUpperBnd)[i] - u[i];
				for(int cmp = 0; cmp < blockSize; cmp++)
				{
					if(((flags >> cmp) & 1)
						&& BlockRef(t, cmp) < BlockRef(diff, cmp))
					{
						// set the "Dirichlet row" and recompute the correction
						for(size_t k = 0; k < blockSize; k++) BlockRef(A, i, k) = 0;
						BlockRef(A, i, i) = 1;
						BlockRef(s, i) = BlockRef(diff, i);
						InverseMatMult(t, 1.0, A, s);
					}
				}
			}
		}
		
		c[i] = t;
	}
	
	return true;
}

/**
 * The general iteration loop of the "non-linear linear solver"
 */
template <typename TAlgebra>
bool ObstacleGaussSeidel<TAlgebra>::apply
(
	vector_type& u, ///< the solution to compute
	const vector_type& b ///< the right-hand side
)
{
	#ifdef UG_PARALLEL
	if(!b.has_storage_type(PST_ADDITIVE) || !u.has_storage_type(PST_CONSISTENT))
		UG_THROW("ObstacleGaussSeidel::apply: Inadequate parallel storage format of Vectors: "
					<< b.get_storage_type() << " for b (expected " << PST_ADDITIVE << "), "
					<< u.get_storage_type() << " for u (expected " << PST_CONSISTENT << ")");
	#endif

//	debug output
	if(this->vector_debug_writer_valid())
		write_debug(b, std::string("CGS_RHS") + ".vec");
	
// 	create correction
	SmartPtr<vector_type> spC = u.clone_without_values();
	vector_type& c = *spC;
	#ifdef UG_PARALLEL
		// this is ok if clone_without_values() inits with zeros
		c.set_storage_type(PST_CONSISTENT);
	#endif
	
// 	create defect
	SmartPtr<vector_type> spD = b.clone();
	vector_type& d = *spD;
	#ifdef UG_PARALLEL
		// this is ok if clone_without_values() inits with zeros
		d.set_storage_type(PST_ADDITIVE);
	#endif
	
	int loopCnt = 0;
	
// 	build defect:  d := b - A*u
	linear_operator()->apply_sub(d, u); // initially, d = b
	
//	compute the first correction
	enter_precond_debug_section(loopCnt);
	if(! compute_correction(c, d, u))
	{
		this->leave_vector_debug_writer_section();
		return false;
	}
	this->leave_vector_debug_writer_section();

	prepare_conv_check();
	convergence_check()->start(c); // we use correction instead of the defect in the termination criterion

	write_debugSolCorDef(u, c, d, loopCnt, false);

// 	Iteration loop
	while(!convergence_check()->iteration_ended())
	{
	// 	add correction to solution
		u += c;
		
	// 	build defect:  d := b - A*u
		d = b;
		linear_operator()->apply_sub(d, u);
		write_debugSolCorDef(u, c, d, ++loopCnt, true);
	
	//	compute the new correction
		enter_precond_debug_section(loopCnt);
		if(! compute_correction(c, d, u))
		{
			this->leave_vector_debug_writer_section();
			return false;
		}
		this->leave_vector_debug_writer_section();

	// 	compute norm of new defect (in parallel)
		convergence_check()->update(c);
	}

//	write some information when ending the iteration
	if(!convergence_check()->post())
	{
		UG_LOG("ERROR in 'ObstacleGaussSeidel::apply': post-convergence-check signaled failure. Aborting.\n");
		return false;
	}

//	we're done
	return true;
}

/**
 * Applies the method and converts the right-hand side b to the defect
 */
template <typename TAlgebra>
bool ObstacleGaussSeidel<TAlgebra>::apply_return_defect
(
	vector_type& u, ///< solution to get
	vector_type& d ///< initially the rhs, finally - the defect
)
{
//	solve u
	if(!apply(u, d)) return false;

//	update defect
	linear_operator()->apply_sub(d, u);

//	we're done
	return true;
}

} //end of namespace SmallStrainMechanics
} //end of namespace ug

/* End of File */