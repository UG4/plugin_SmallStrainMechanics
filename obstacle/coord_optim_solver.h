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

#ifndef __SMALL_STRAIN_MECH_COORD_OPTIM_H__
#define __SMALL_STRAIN_MECH_COORD_OPTIM_H__
#include <iostream>
#include <string>

#include "lib_algebra/operator/interface/matrix_operator_inverse.h"
#include "lib_algebra/operator/interface/linear_solver_profiling.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{
namespace SmallStrainMechanics{

/// class of the quadratic coordinate-direction optimizer (a type of the Gauss-Seidel method)
/**
 * The method is based on the optimization along the coordinate directions
 * (a Gauss-Seidel-like method) with constraints (obstacles).
 *
 * \tparam 		TAlgebra		algebra type
 */
template <typename TAlgebra>
class ObstacleGaussSeidel
	: public IMatrixOperatorInverse<typename TAlgebra::matrix_type, typename TAlgebra::vector_type>,
	  public DebugWritingObject<TAlgebra>
{
	typedef IMatrixOperatorInverse<typename TAlgebra::matrix_type, typename TAlgebra::vector_type> base_type;

protected:
	using base_type::convergence_check;
	using base_type::linear_operator;
	
	using DebugWritingObject<TAlgebra>::debug_writer;
	using DebugWritingObject<TAlgebra>::write_debug;

public:
///	Algebra type
	typedef TAlgebra algebra_type;
	static const int blockSize = algebra_type::blockSize;

///	Vector type
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename vector_type::flag_unit_type flag_unit_type;

///	Matrix type
	typedef typename TAlgebra::matrix_type matrix_type;
	
///	Flag manager type
	typedef DoFFlagManager<algebra_type> flag_mng_type;
	static const flag_unit_type flagBlockFrame = flag_mng_type::flagBlockFrame;
	
public:

///	constructors
	ObstacleGaussSeidel() {}

public:
///	returns the name of the solver
	virtual const char* name() const {return "Coordinate Direction Optimizer";}

///	returns if parallel solving is supported
	virtual bool supports_parallel() const {return true;}
	
///	initializer
	virtual bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
	{
		base_type::m_spLinearOperator = Op;
		m_pA = & Op->get_matrix ();
		return true;
	}

///	solves the system and returns the last defect
	virtual bool apply(vector_type& u, const vector_type& b);
	
/// applies the method and converts the right-hand side b to the defect if needed
	virtual bool apply_return_defect(vector_type& u, vector_type& b);

///	parses the contact bnd flag
	void set_contact_bnd_flag
	(
		const char* bnd_flag ///< names of the flag
	)
	{
		if(! flag_mng_type::get()->get_flag_offsets(bnd_flag, m_obst_bnd_flag[0], m_obst_bnd_flag[1]))
			UG_THROW("ObstacleGaussSeidel: Flag '" << bnd_flag << "' not found");
		m_obst_bnd_flag_name = bnd_flag;
	}
	
///	sets the lower bound
	void set_lower_bnd
	(
		const vector_type * lower ///< the vector with the lower constraints
	)
	{
		if (lower == NULL)
			m_spLowerBnd = SPNULL;
		else
			m_spLowerBnd = lower->clone();
	}
	
///	sets the upper bound
	void set_upper_bnd
	(
		const vector_type * upper ///< the vector with the lower constraints
	)
	{
		if (upper == NULL)
			m_spUpperBnd = SPNULL;
		else
			m_spUpperBnd = upper->clone();
	}
	
protected:

///	computes one step of the constrained Gauss-Seidel iteration
	bool compute_correction(vector_type &c, const vector_type &d, const vector_type &u);

///	prepares the convergence check output
	void prepare_conv_check()
	{
		convergence_check()->set_name(name());
		convergence_check()->set_symbol('%');
	}

/// debugger output: solution, correction, defect
	void write_debugSolCorDef(vector_type &x, vector_type &c, vector_type &d, int loopCnt, bool bWriteC)
	{
		if(!this->vector_debug_writer_valid()) return;
		char ext[20]; snprintf(ext, sizeof(ext),"_iter%03d", loopCnt);
		write_debug(d, std::string("CGS_Defect") + ext + ".vec");
		if(bWriteC) write_debug(c, std::string("LS_Correction") + ext + ".vec");
		write_debug(x, std::string("CGS_Solution") + ext + ".vec");
	}
	
/// debugger section for the preconditioner
	void enter_precond_debug_section(int loopCnt)
	{
		if(!this->vector_debug_writer_valid()) return;
		char ext[20]; snprintf(ext, sizeof(ext),"_iter%03d", loopCnt);
		this->enter_vector_debug_writer_section(std::string("CGS_Precond_") + ext);
	}
	
private:

	const matrix_type* m_pA; ///< the system matrix
	
	std::string m_obst_bnd_flag_name; ///< the name of the constraint flag
	size_t m_obst_bnd_flag[2];	///< constraint flag indices
	SmartPtr<vector_type> m_spLowerBnd; ///< constraints from below (a cloned vector)
	SmartPtr<vector_type> m_spUpperBnd; ///< constraints from above (a cloned vector)
};

} //end of namespace SmallStrainMechanics
} // end namespace ug

#include "coord_optim_solver_impl.h"

#endif /* __SMALL_STRAIN_MECH_COORD_OPTIM_H__ */

/* End of File */

