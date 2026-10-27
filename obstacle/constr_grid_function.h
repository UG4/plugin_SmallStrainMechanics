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

// ug4 headers
#include "lib_disc/function_spaces/grid_function.h"

#ifndef __SMALL_STRAIN_MECH_CONSTR_GRID_FUNCTION_H__
#define __SMALL_STRAIN_MECH_CONSTR_GRID_FUNCTION_H__

namespace ug{
namespace SmallStrainMechanics{

/// class for grid functions with upper and lower constraints
/**
 * This class adds two optional vectors to a grid function: the lower and the
 * upper bounds.
 *
 * IMPORTANT REMARK:
 * This implementation does not include the adaption of the bounds.
 * It should be extended with the "managing" facilities.
 *
 * \tparam 	TDomain				domain type
 * \tparam	TAlgebra			algebra type
 */
template <typename TDomain, typename TAlgebra>
class ConstrainedGridFunction
:	public GridFunction<TDomain, TAlgebra>
{
public:
	typedef TDomain domain_type; ///< domain type
	typedef TAlgebra algebra_type; ///< algebra type
	
	typedef typename algebra_type::vector_type vector_type; //< vector type
	
	typedef ConstrainedGridFunction<TDomain, TAlgebra> this_type;
	typedef GridFunction<TDomain, TAlgebra> base_type;
	
public:

///	Inherit constructors from GridFunction
    using base_type::base_type;

///	Copy constructor: copies the values and the bounds (if any)
    ConstrainedGridFunction(const this_type& v) : base_type(v)
    {
        if(v.m_spLowerBnd.valid())
        	m_spLowerBnd = v.m_spLowerBnd->clone();
        if(v.m_spUpperBnd.valid())
        	m_spLowerBnd = v.m_spUpperBnd->clone();
    }

/// virtual clone using covariant return type
	virtual this_type* virtual_clone() const {return new this_type(*this);}
	
/// virtual clone using covariant return type excluding values
	virtual this_type* virtual_clone_without_values() const
	{
		this_type* p = new this_type(this->m_spApproxSpace, this->m_spDD, this->m_bManaged);
		p->enable_redistribution(this->redistribution_enabled());
		if(p->size() != this->size())
			p->resize(this->size());
#ifdef UG_PARALLEL
		p->set_layouts(this->layouts());
#endif
        if(m_spLowerBnd.valid())
        	p->m_spLowerBnd = m_spLowerBnd->clone();
        if(m_spUpperBnd.valid())
        	p->m_spLowerBnd = m_spUpperBnd->clone();
	
		return p;
	}
	
public:

///	get the vector of the lower bound
	SmartPtr<vector_type> lower_bnd() {return m_spLowerBnd;}

///	get the vector of the upper bound
	SmartPtr<vector_type> upper_bnd() {return m_spUpperBnd;}
	
///	create the lower bound without values
	void create_lower_bnd()
	{
		if(m_spLowerBnd.valid()) return;
		m_spLowerBnd = make_sp(new vector_type(this->size()));
	}

///	create the upper bound without values
	void create_upper_bnd()
	{
		if(m_spUpperBnd.valid()) return;
		m_spUpperBnd = make_sp(new vector_type(this->size()));
	}

public:

	SmartPtr<vector_type> m_spLowerBnd; ///< vector for the lower bound
	SmartPtr<vector_type> m_spUpperBnd; ///< vector for the upper bound
};

} //end of namespace SmallStrainMechanics
} //end of namespace ug

#endif // __SMALL_STRAIN_MECH_CONSTR_GRID_FUNCTION_H__

/* End of File */
