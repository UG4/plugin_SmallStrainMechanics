/*
 * hooke.h
 *
 *  Created on: 03.02.2014
 *      Author: raphaelprohl
 */

#ifndef HOOKE_H_
#define HOOKE_H_

namespace ug{
namespace SmallStrainMechanics{

/// Hookes Law
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
 * <li> D. Braess. Finite Elemente. Theorie, schnelle Löser und Anwendungen in der
 * <li> Elastizitätstheorie, Springer
 * <ul>
 */

template <typename TDomain>
class HookeLaw
	: public IMaterialLaw<TDomain>
{
	private:
	///	Base class type
		typedef IMaterialLaw<TDomain> base_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	constructor
		HookeLaw(): IMaterialLaw<TDomain>()
		{
			base_type::m_bConstElastTens = true;
		};

	///	Destructor
		~HookeLaw(){};

	public:
	////////////////////////////
	// INTERFACE-METHODS
	////////////////////////////
		void init();

	///	computes the cauchy stress tensor sigma at an integration point 'ip'
		void stressTensor(MathMatrix<dim,dim>& stressTens, const size_t ip,
				const MathMatrix<dim, dim>& GradU);

	///	computes the elasticity tensor; commonly denoted by C
		virtual inline SmartPtr<MathTensor4<TDomain::dim,TDomain::dim,TDomain::dim,TDomain::dim> >
				elasticityTensor();

		virtual SmartPtr<MathTensor4<TDomain::dim,TDomain::dim,TDomain::dim,TDomain::dim> >
			elasticityTensor(const size_t ip, const MathMatrix<dim, dim>& GradU){
			return elasticityTensor();};


	public:
	///	set elasticity tensor for orthotropic materials
		void set_elasticity_tensor_orthotropic(
				const number C11, const number C12, const number C13,
							const number C22, const number C23,
										const number C33,
												const number C44,
														const number C55,
																const number C66 );

	///	set hooke elasticity tensor for isotropic materials, (in 2D: plane-strain-case)
		void set_hooke_elasticity_tensor(const number lambda, const number mu);
		void set_hooke_elasticity_tensor_E_nu(const number E, const number nu);

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

#include "hooke_impl.h"

#endif /* HOOKE_H_ */
