/*
 * Copyright (c) 2019: Ruhr University Bochum
 * Author: Andreas Vogel
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

#ifndef SCALED_HOOKE_LAW_H_
#define SCALED_HOOKE_LAW_H_

#include "mat_law_interface.h"
#include "hooke.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_algebra/cpu_algebra_types.h"



namespace ug{
namespace SmallStrainMechanics{

/// \addtogroup small_strain_mechanics
/// \{




template <typename TDomain>
class IScaledHookeLaw : public HookeLaw<TDomain>
{
	private:
	///	Base class type
		typedef HookeLaw<TDomain> base_type;

	///	own type
		typedef IScaledHookeLaw<TDomain> this_type;


	public:
	///	World dimension
		static const int dim = base_type::dim;

	///	base element type
		typedef typename base_type::TBaseElem TBaseElem;

	public:
	///	constructor
		IScaledHookeLaw(	SmartPtr< GridFunction<TDomain, CPUAlgebra> > spScaling,
							SmartPtr< GridFunction<TDomain, CPUAlgebra> > spEnergy)
			: HookeLaw<TDomain>(), 
			m_spScaling(spScaling), m_spEnergy(spEnergy), 
			m_pScaling_elem(NULL), m_pEnergy_elem(NULL) 
		{}

	///	Destructor
		~IScaledHookeLaw(){};

	/// reset values explicitly
		void init_internal_vars(TBaseElem* elem, const size_t numIP)
		{
			if (!base_type::m_bInit){
				base_type::m_bInit = true;

				m_spScaling->set(1.0);
				m_spEnergy->set(0.0);
			}
		}

		void internal_vars(TBaseElem* elem)
		{
			const size_t fct = 0; // \todo: generalize
			std::vector<DoFIndex> ind;
			if(m_spScaling->inner_dof_indices(elem, fct, ind) != 1)
				UG_THROW("Wrong number dofs");

			m_pScaling_elem = & DoFRef(*m_spScaling, ind[0]);
			m_pEnergy_elem = &DoFRef(*m_spEnergy, ind[0]);
		}

		virtual void post_process_energy_on_curr_elem(){}

	public:
		virtual number scaling_on_curr_elem() {return *m_pScaling_elem;}
		number& energy_on_curr_elem() {return *m_pEnergy_elem;}

	protected:
		SmartPtr<GridFunction<TDomain, CPUAlgebra > > m_spScaling;
		SmartPtr<GridFunction<TDomain, CPUAlgebra > > m_spEnergy;

		number* m_pScaling_elem;
		number* m_pEnergy_elem;
};




/// Material Law: 
/**
 * 	
 */
template <typename TDomain>
class DamageLaw
	: public IScaledHookeLaw<TDomain>
{
	private:
	///	Base class type
		typedef IScaledHookeLaw<TDomain> base_type;

	///	own type
		typedef DamageLaw<TDomain> this_type;

	protected:
		using base_type::m_pEnergy_elem;
		using base_type::m_pScaling_elem;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	///	base element type
		typedef typename base_type::TBaseElem TBaseElem;

	public:
	///	constructor
		DamageLaw(	SmartPtr< GridFunction<TDomain, CPUAlgebra> > spF,
					SmartPtr< GridFunction<TDomain, CPUAlgebra> > spPsi0)
			: IScaledHookeLaw<TDomain>(spF, spPsi0)
		{}

	///	Destructor
		~DamageLaw() {};


		virtual void post_process_energy_on_curr_elem()
		{
			(*m_pEnergy_elem) /= (*m_pScaling_elem);
		}
};

/// Material Law: 
/**
 * 	
 */
template <typename TDomain>
class TopologyOptimLaw
	: public IScaledHookeLaw<TDomain>
{
	private:
	///	Base class type
		typedef IScaledHookeLaw<TDomain> base_type;

	///	own type
		typedef TopologyOptimLaw<TDomain> this_type;

	protected:
		using base_type::m_pEnergy_elem;
		using base_type::m_pScaling_elem;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	///	base element type
		typedef typename base_type::TBaseElem TBaseElem;

	public:
	///	constructor
		TopologyOptimLaw(	SmartPtr< GridFunction<TDomain, CPUAlgebra> > spChi,
							SmartPtr< GridFunction<TDomain, CPUAlgebra> > spDrivingForce,
							int expPenalize)
			: IScaledHookeLaw<TDomain>(spChi, spDrivingForce),
			m_expPenalize(expPenalize)
		{}

	///	Destructor
		~TopologyOptimLaw() {};

		virtual number scaling_on_curr_elem() {return std::pow(*m_pScaling_elem, m_expPenalize);}


		virtual void post_process_energy_on_curr_elem()
		{
			(*m_pEnergy_elem) = ((m_expPenalize-1) * (*m_pEnergy_elem)) / (*m_pScaling_elem);
		}

	protected:
		int m_expPenalize;
};

}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#endif /* SCALED_HOOKE_LAW_H_ */
