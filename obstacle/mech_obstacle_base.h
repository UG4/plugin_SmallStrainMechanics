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

/**
 * Interfaces for the mechanical obstacle conditions for the small strain mechanics
 * discretizations.
 */

#ifndef __SMALL_STRAIN_MECH_MECH_OBSTACLE_H__
#define __SMALL_STRAIN_MECH_MECH_OBSTACLE_H__

#include <vector>
#include <map>

// ug4 hearders
#include "lib_algebra/dof_flag_manager.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"

namespace ug{
namespace SmallStrainMechanics{

/// Class for local basis based constraints
/**
 * This is a base class for constraints based on the local basis,
 * e.g. the normal and the tangential directions
 */
template <typename TDomain>
class SignoriniObstacle
{
public:

/// own type
	typedef SignoriniObstacle<TDomain> this_type;
	
///	domain type
	typedef TDomain domain_type;

///	world dimension
	static const int dim = domain_type::dim;
	
///	constraint value type
	typedef UserData<number, dim> constr_value_type;
	
///	constraint direction type
	typedef UserData<MathVector<dim>, dim> constr_dir_type;

public:
	
///	Dummy constructor
	SignoriniObstacle() : m_value(NULL), m_normal(NULL) {}
	
///	Constructor
	SignoriniObstacle
	(
		SmartPtr<constr_value_type> value,
		SmartPtr<constr_dir_type> normal,
		std::string name
	)
	: m_value(value), m_normal(normal), m_name(name)
	{}
	
///	Returns the UserData of the normal
	SmartPtr<constr_dir_type> normal() {return m_normal;}
	
///	Returns the UserData of the constraint value
	SmartPtr<constr_value_type> value() {return m_value;}
	
private:
	
	SmartPtr<constr_value_type> m_value; ///< constraint value
	SmartPtr<constr_dir_type> m_normal; ///< constraint normal direction

///	Name of the constraint
	std::string m_name;
};

/// Class of the constraints manager
/**
 * This class collects all the constraints and provides the interface for
 * the discretization.
 *
 * \tparam TDomain	domain type
 * \tparam TAlgebra	algebra type
 */
template <typename TDomain, typename TAlgebra>
class SignoriniConstraint
	: public IDomainConstraint<TDomain, TAlgebra>
{
public:
/// own type
	typedef SignoriniConstraint<TDomain, TAlgebra> this_type;
	
///	domain type
	typedef TDomain domain_type;

///	world dimension
	static const int dim = domain_type::dim;
	
///	type of the position vectors
	typedef typename domain_type::position_type position_type;

///	type of algebra
	typedef TAlgebra algebra_type;

///	type of algebra matrix
	typedef typename algebra_type::matrix_type matrix_type;

///	type of algebra vector
	typedef typename algebra_type::vector_type vector_type;
	
///	the type map for the constrants
	typedef typename std::map<int, SmartPtr<SignoriniObstacle<domain_type> > > constr_map_type;

public:

///	Constructor
	SignoriniConstraint
	(
		SmartPtr<domain_type> domain, ///< the domain
		const char* fct_names, ///< names of the functions
		const char* flag_name ///< name of the flag to mark the boundary
	)
	:	m_spDomain(domain)
	{
		set_fcts_names(fct_names);
		set_flag(flag_name);
	}
	
///	returns the associated domain
	SmartPtr<domain_type> domain() const {return m_spDomain;}
	
///	remove all the constraints
	void clear () {m_mspConstraints.clear(); m_vsFctNames.clear();}
	
///	set function names
	void set_fcts_names
	(
		std::string fctNames ///< (comma separated) function names
	)
	{
		TokenizeString(fctNames, m_vsFctNames);
		for(size_t i = 0; i < m_vsFctNames.size(); i++)
			RemoveWhitespaceFromString(m_vsFctNames[i]);
	}
	
///	parses the flags
	void set_flag
	(
		std::string bnd_flag ///< names of the flag
	)
	{
		if(! DoFFlagManager<algebra_type>::get()->get_flag_offsets(bnd_flag.c_str(),
			m_bnd_flag[0], m_bnd_flag[1]))
			UG_THROW("SignoriniConstraint: Flag '" << bnd_flag << "' not found");
	}
	
///	add a new constraint for a subset index
	void add
	(
		int ssi, ///< subset index
		SmartPtr<SignoriniObstacle<domain_type> > spConstr ///< the constraint
	)
	{
		if(m_mspConstraints.find(ssi) != m_mspConstraints.end())
			UG_THROW("SignoriniConstraint: Attempt to add the second constraint to subset id " << ssi);
		m_mspConstraints[ssi] = spConstr;
	}
	
///	add a new constraint by subset names
	void add
	(
		std::string ss_names, ///< list of subset names
		SmartPtr<SignoriniObstacle<domain_type> > spConstr ///< the constraint
	)
	{
		std::vector<std::string> vNames;
		SubsetGroup ssg;
		
		TokenizeString(ss_names, vNames);
		for(size_t i = 0; i < vNames.size(); i++)
			RemoveWhitespaceFromString(vNames[i]);
		ssg.set_subset_handler(m_spDomain->subset_handler());
		ssg.add(vNames);
		
		for(size_t i = 0; i < ssg.size(); i++)
			add(ssg[i], spConstr);
	}
	
//	IDomainConstraint interface
	
///	adapts jacobian to enforce constraints
	virtual void adjust_jacobian(matrix_type& J, const vector_type& u,
								 ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
								 ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
								 const number s_a0 = 1.0)
	{
		matrix_type transform;
		vector_type bnd_gf;
		
		transform.resize_and_clear(J.num_rows(), J.num_cols());
		bnd_gf.resize(J.num_cols());
		
		FunctionGroup fctGrp (dd->function_pattern(), m_vsFctNames);
		
		mark_signorini_bnd(*dd, fctGrp, transform, bnd_gf, time);
		transform_matrix(transform, J);
	}

///	adapts defect to enforce constraints
	virtual void adjust_defect(vector_type& d, const vector_type& u,
							   ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
							   ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
							   const std::vector<number>* vScaleMass = NULL,
							   const std::vector<number>* vScaleStiff = NULL)
	{
		matrix_type transform;
		vector_type bnd_gf;
		
		transform.resize_and_clear(d.size(), d.size());
		bnd_gf.resize(d.size());
		
		FunctionGroup fctGrp (dd->function_pattern(), m_vsFctNames);
		
		mark_signorini_bnd(*dd, fctGrp, transform, bnd_gf, time);
		transform_vector(transform, d);
	}

///	adapts matrix and rhs (linear case) to enforce constraints
	virtual void adjust_linear(matrix_type& mat, vector_type& rhs,
							   ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0)
	{
		matrix_type transform;
		vector_type bnd_gf;
		
		transform.resize_and_clear(mat.num_rows(), mat.num_cols());
		bnd_gf.resize(mat.num_cols());
		
		FunctionGroup fctGrp (dd->function_pattern(), m_vsFctNames);
		
		mark_signorini_bnd(*dd, fctGrp, transform, bnd_gf, time);
		transform_matrix(transform, mat);
		transform_vector(transform, rhs);
	}

///	adapts a rhs to enforce constraints
	virtual void adjust_rhs(vector_type& rhs, const vector_type& u,
							ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0)
	{
		matrix_type transform;
		vector_type bnd_gf;
		
		transform.resize_and_clear(rhs.size(), rhs.size());
		bnd_gf.resize(rhs.size());
		
		FunctionGroup fctGrp (dd->function_pattern(), m_vsFctNames);
		
		mark_signorini_bnd(*dd, fctGrp, transform, bnd_gf, time);
		transform_vector(transform, rhs);
	}

///	sets the constraints in a solution vector
	virtual void adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd, int type,
								 number time = 0.0)
	{
		matrix_type transform;
		vector_type bnd_gf;
		
		transform.resize_and_clear(u.size(), u.size());
		bnd_gf.resize(u.size());
		
		FunctionGroup fctGrp (dd->function_pattern(), m_vsFctNames);
		
		mark_signorini_bnd(*dd, fctGrp, transform, bnd_gf, time);
		transform_vector(transform, u);
	}

///	returns the type of the constraints
	virtual int type() const {return CT_ASSEMBLED;} //ToDo: Is it a correct type here?

protected:

///	returns the constraint map
	constr_map_type& constr_maps() {return m_mspConstraints;}
	
//	The assembling interface

private:

/// marks the Signorini boundary DoFs and computes the normal vectors for a particular element type
	template <typename TBaseElem> 
	void mark_signorini_bnd_elem
	(
		const DoFDistribution& dd, ///< DoF distribution of the solution
		FunctionGroup& fctGrp, ///< function group of the displacement
		matrix_type& transform, ///< the transformation (rotation) matrix
		vector_type& bnd_gf, ///< distribution of the bound
		number time ///< the time argument
	);

protected:
	
/// marks the Signorini boundary DoFs and computes the normal vectors
	void mark_signorini_bnd
	(
		const DoFDistribution& dd, ///< DoF distribution of the solution
		FunctionGroup& fctGrp, ///< function group of the displacement
		matrix_type& transform, ///< the transformation (rotation) matrix
		vector_type& bnd_gf, ///< distribution of the bound
		number time ///< the time argument
	);
	
///	apply the transformation to a matrix
	void transform_matrix
	(
		matrix_type& transform, ///< the transformation (rotation) matrix
		matrix_type& A ///< the matrix to transform
	)
	{
		matrix_type A_t;
		CreateAsMultiplyOf(A_t, transform, A, transform);
		A = A_t;
	}
	
///	apply the transformation to a vector
	void transform_vector
	(
		matrix_type& transform, ///< the transformation (rotation) matrix
		vector_type& v ///< the vector to transform
	)
	{
		vector_type v_t;
		v_t.resize(v.size());
		transform.axpy(v_t, 0.0, v_t, 1.0, v);
	}
	
private:
	
	SmartPtr<domain_type> m_spDomain;
	std::vector<std::string> m_vsFctNames;
	constr_map_type m_mspConstraints;

	size_t m_bnd_flag[2];	///< flag indices
};

} //end of namespace SmallStrainMechanics
} //end of namespace ug

// Implementation of the functions
#include "mech_obstacle_base_impl.h"

#endif // __SMALL_STRAIN_MECH_MECH_OBSTACLE_H__

/* End of File */
