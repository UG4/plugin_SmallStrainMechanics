/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef SKIN_LAW_H_
#define SKIN_LAW_H_

namespace ug{
namespace SmallStrainMechanics{

/// \addtogroup small_strain_mechanics
/// \{
/// Material Law: Hookes Law modelling linear elastic behavior for isotropic, homogenous materials
/**
 * 	The classical material law to describe linear elastic behavior for
 * 	isotropic, homogeneous materials. This law is based on the assumption
 * 	of small strains, since only the linearized strain tensor
 * 	(commonly denoted with \epsilon) is used. With $u$ being the
 * 	displacement vector, the strain tensor \epsilon is defined by
 *
 * 		\epsilon = \frac{1}{2}(\nabla u + (\nabla u)^T),
 *
 * 	i.e. higher orders of the gradient of u are neglected here
 * 	(-> 'linearized strain tensor').
 *
 *	Additionally Hooke`s law assumes a linear connection of stresses
 *	and strains. Therefore, the stresses \sigma are described by a
 *	law of the form,
 *
 *		\sigma = C : \epsilon
 *
 *	with C being a constant tensor of rank 4, called 'ElasticityTensor'.
 *	\sigma denotes the cauchy-stress tensor, which is, just like the
 *	strain tensor \epsilon, a symmetric tensor of rank 2.
 *	':' denotes the tensor contraction, in components,
 *
 *	\sigma_{ij} = C_{ijkl} \epsilon{kl}.
 *
 * \tparam 	int		dim
 *
 * References:
 * <ul>
 * <li> D. Braess. Finite Elemente. Theorie, schnelle L�ser und Anwendungen in der
 * <li> Elastizit�tstheorie, Springer
 * <ul>
 */

template <typename TDomain>
class SkinMaterialLaw
	: public IMaterialLaw<TDomain>
{
	private:
	///	Base class type
		typedef IMaterialLaw<TDomain> base_type;

	public:
	///	World dimension
		static const int dim = TDomain::dim;

	public:
	///	constructor
		SkinMaterialLaw(): IMaterialLaw<TDomain>()
		{
			base_type::m_bConstElastTens = false;
		};

	///	Destructor
		~SkinMaterialLaw(){};

	public:
	////////////////////////////
	// INTERFACE-METHODS
	////////////////////////////
		void init();

	///	computes the cauchy stress tensor sigma at an integration point 'ip'
		void stressTensor(MathMatrix<dim,dim>& stressTens, const size_t ip,
							const MathMatrix<dim, dim>& GradU){
			UG_THROW("SkinMaterialLaw: only position dependent material law supported");
		}
		void stressTensor(MathMatrix<dim,dim>& stressTens, const size_t ip,
					const MathVector<dim>& x, const MathMatrix<dim, dim>& GradU);

	///	computes the elasticity tensor; commonly denoted by C
		virtual inline SmartPtr<MathTensor4<dim,dim,dim,dim> > elasticityTensor(){
			UG_THROW("SkinMaterialLaw: only position dependent material law supported");
		}

		virtual SmartPtr<MathTensor4<dim,dim,dim,dim> >
			elasticityTensor(const size_t ip, const MathMatrix<dim, dim>& GradU){
			UG_THROW("SkinMaterialLaw: only position dependent material law supported");
		}

		virtual SmartPtr<MathTensor4<TDomain::dim,TDomain::dim,TDomain::dim,TDomain::dim> >
			elasticityTensor(const size_t ip, const MathVector<dim>& x, const MathMatrix<dim, dim>& GradU);


	public:
	///	set hooke elasticity tensor for isotropic materials, (in 2D: plane-strain-case)
	/// \{
		void set_hooke_elasticity_tensor(const number lambda, const number mu);
		void set_hooke_elasticity_tensor_E_nu(const number E, const number nu);

		void set_hooke_elasticity_tensor_plain_stress_E_nu(const number E, const number nu);
		void set_hooke_elasticity_tensor_plain_strain_E_nu(const number E, const number nu);
	/// \}

	public:
		using base_type::m_materialConfiguration;

	private:
		void strainTensor(MathMatrix<dim,dim>& strainTens, const MathMatrix<dim, dim>& GradU);

	private:
	/// elasticity tensor
		SmartPtr<MathTensor4<dim,dim,dim,dim> > m_spElastTensorFunct;


};

}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#include "skin_law_impl.h"

#endif /* SKIN_LAW_H_ */
