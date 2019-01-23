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

#ifndef DAMAGE_LAW_H_
#define DAMAGE_LAW_H_

#include "mat_law_interface.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_algebra/cpu_algebra_types.h"



namespace ug{
namespace SmallStrainMechanics{

/// \addtogroup small_strain_mechanics
/// \{
/// Material Law: Prandtl-Reuss law modelling elastoplastic material behavior where the elastic
///	part is considered as linear. The flow-condition is of von-Mises-type.
/**
 * 	This class implements a material law for small strain elastoplastic material behavior
 *
 * 	It is supposed, that the linearized strain tensor could be decomposed additively:
 *
 *  eps = eps_e + eps_p.
 *
 *  The plastic behavior is described by a flow-condition and a flow-rule for the plastic
 *  evolution (\frac{\partial eps_p){\partial t} = ...). The flow-condition is of
 *  von-Mises-type and the flow-rule is associative. To treat the plastic equations
 *  we use the well-established return-mapping-algorithm. Its classical form is valid for the
 *  3d-case and the plane strain-case, but not for the plane stress-case!
 *
 * References:
 * <ul>
 * <li> J.C. Simo and T.J.R. Hughes. Computational Inelasticity. Springer, New York (1998), chapter 3.3.1
 * <li>
 * <li> F.-J. Barthold, M. Schmidt and E. Stein. Error indicators and mesh refinements for
 * <li> finite-element computations of elastoplastic deformations. Computational Mechanics Vol. 22, 225-238 (1998)
 *</ul>
 *
 * \tparam TDomain
 */

template <typename TDomain>
class DamageLaw
	: public IMaterialLaw<TDomain>
{
	private:
	///	Base class type
		typedef IMaterialLaw<TDomain> base_type;

	///	own type
		typedef DamageLaw<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	///	base element type
		typedef typename base_type::TBaseElem TBaseElem;

	public:
	///	constructor
//		DamageLaw();
		DamageLaw(	SmartPtr< GridFunction<TDomain, CPUAlgebra> > spPsi0,
					SmartPtr< GridFunction<TDomain, CPUAlgebra> > f);

	///	Destructor
		~DamageLaw(){};

	public:
	////////////////////////////
	// INTERFACE-METHODS
	////////////////////////////
		void init();

	///	computes the cauchy stress tensor sigma at an integration point 'ip'
		void stressTensor(MathMatrix<dim,dim>& stressTens, const size_t ip,
				const MathMatrix<dim, dim>& GradU);

	///	computes the elasticity tensor; commonly denoted by C
		SmartPtr<MathTensor4<TDomain::dim,TDomain::dim,TDomain::dim,TDomain::dim> >
			elasticityTensor(const size_t ip, const MathMatrix<dim, dim>& GradU);

		virtual bool needs_to_add_jac_m(){return false;}


		virtual void init_internal_vars(TBaseElem* elem, const size_t numIP);
		virtual void internal_vars(TBaseElem* elem);
		virtual void update_internal_vars(const size_t ip, const MathMatrix<dim, dim>& GradU){};
		virtual void attach_internal_vars(typename TDomain::grid_type& grid);
		virtual void clear_attachments(typename TDomain::grid_type& grid);

		number& f_on_curr_elem() {return *m_pF_elem;}
		number& psi0_on_curr_elem() {return *m_pPsi0_elem;}

	public:
	///	set hooke elasticity tensor for isotropic materials, (in 2D: plane-strain-case)
	/// \{
		void set_hooke_elasticity_tensor(const number lambda, const number mu);
		void set_hooke_elasticity_tensor_E_nu(const number E, const number nu);

		void set_hooke_elasticity_tensor_plain_stress_E_nu(const number E, const number nu);
		void set_hooke_elasticity_tensor_plain_strain_E_nu(const number E, const number nu);
	/// \}

	///	set elasticity tensor for orthotropic materials
	/// \{
		void set_elasticity_tensor_orthotropic(
			const number C11, 	const number C12, 	const number C13,
								const number C22, 	const number C23,
													const number C33,
																	const number C44,
																		const number C55,
																			const number C66 );

		void set_elasticity_tensor_orthotropic_E_G_nu(
				const number E1, const number E2, const number E3, 
				const number G12, const number G13, const number G23, 
				const number v12, const number v13, const number v23);

		void set_elasticity_tensor_orthotropic_plain_stress_E_G_nu(
			const number E1, const number E2, 
			const number G12, 
			const number v12);

		void set_elasticity_tensor_orthotropic_plain_strain_E_G_nu(
				const number E1, const number E2, const number E3, 
				const number G12, const number G13, const number G23, 
				const number v12, const number v13, const number v23
				);
	/// \}

	public:
		using base_type::m_materialConfiguration;

	public:
		void strainTensor(MathMatrix<dim,dim>& strainTens, const MathMatrix<dim, dim>& GradU);

	private:
	/// elasticity tensor
		SmartPtr<MathTensor4<dim,dim,dim,dim> > m_spElastTensorFunct;


	private:
		SmartPtr<GridFunction<TDomain, CPUAlgebra > > m_spPsi0;
		SmartPtr<GridFunction<TDomain, CPUAlgebra > > m_spF;

		number* m_pF_elem;
		number* m_pPsi0_elem;

/*
	//	std-vector of InternalVars
		struct ElemData{
			number f;
			number psi0;
		};

		ElemData* m_pElemData;

	//	attachment type: attachment of ElemData
		typedef Attachment<ElemData> AElemData;
		AElemData m_aElemData;

	//	the attachment accessor
		typedef Grid::AttachmentAccessor<TBaseElem, AElemData>	ElemDataAccessor;
		ElemDataAccessor m_aaElemData;
*/
};

}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#include "damage_law_impl.h"

#endif /* DAMAGE_LAW_H_ */
