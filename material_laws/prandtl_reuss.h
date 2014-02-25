/*
 * prandtl_reuss.h
 *
 *  Created on: 05.02.2014
 *      Author: raphaelprohl
 */

#ifndef PRANDTL_REUSS_H_
#define PRANDTL_REUSS_H_

namespace ug{
namespace SmallStrainMechanics{

/// Prandtl-Reuss Law for ElastoPlasticity
/**
 * 	A material law for small strain elastoplastic material behavior
 *
 * 	It is supposed, that the linearized strain tensor could be decomposed additively:
 *
 *  eps = eps_e + eps_p.
 *
 *  The plastic behavior is described by a flow-condition and a flow-rule for the plastic
 *  evolution (\frac{\partial eps_p){\partial t} = ...). The flow-condition is of
 *  von-Mises-type and the flow-rule is associative. To treat the plasticity
 *  we use the well-established return-mapping-algorithm. Its classical form is valid for the
 *  3d-case and the plane strain-case, but not for the plane stress-case!
 *
 * References:
 * <ul>
 * <li> J. C. Simo and T.J.R. Hughes. Computational Inelasticity. Spinger, New York (1998), chapter 3.3.1
 * <li> F.-J. Barthold, M. Schmidt and E. Stein. Error indicators and mesh refinements for
 * <li> finite-element computations of elastoplastic deformations. Computational Mechanics Vol. 22, 225-238 (1998)
 *</ul>
 *
 * \tparam TDomain
 */

template <typename TDomain>
class PrandtlReuss
	: public IMaterialLaw<TDomain>
{
	private:
	///	Base class type
		typedef IMaterialLaw<TDomain> base_type;

	///	own type
		typedef PrandtlReuss<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	///	base element type
		typedef typename base_type::TBaseElem TBaseElem;

	public:
	///	constructor
		PrandtlReuss();

	///	Destructor
		~PrandtlReuss(){};

	///	set-methods for material constants
		void set_bulk_modulus(const number bulkModulus)
		{matConsts.kappa = bulkModulus;
		std::stringstream ss; ss << m_materialConfiguration << " bulk modulus kappa = " << bulkModulus << "\n";
		m_materialConfiguration = ss.str();}

		void set_shear_modulus(const number shearModulus)
		{matConsts.mu = shearModulus;
		std::stringstream ss; ss << m_materialConfiguration << " shear modulus mu = " << shearModulus << "\n";
		m_materialConfiguration = ss.str();}

		void set_initial_flow_stress(const number initialFlowStress)
		{matConsts.K_0 = initialFlowStress;
		std::stringstream ss; ss << m_materialConfiguration << " initial flow-stress K_0 = " << initialFlowStress << "\n";
		m_materialConfiguration = ss.str();}

		void set_residual_flow_stress(const number resFlowStress)
		{matConsts.K_inf = resFlowStress;
		std::stringstream ss; ss << m_materialConfiguration << " residual flow-stress (saturation stress) K_inf = " << resFlowStress << "\n";
		m_materialConfiguration = ss.str();}

		void set_hardening_modulus(const number hardModulus)
		{matConsts.Hard = hardModulus; m_bHardModulus = true;
		std::stringstream ss; ss << m_materialConfiguration << " hardening modulus = " << hardModulus << "\n";
		m_materialConfiguration = ss.str();}

		void set_hardening_exponent(const number hardExponent)
		{matConsts.omega = hardExponent; m_bHardExp = true;
		std::stringstream ss; ss << m_materialConfiguration << " hardening exponent = " << hardExponent << "\n";
		m_materialConfiguration = ss.str();}

		void set_hardening_behavior(int hard);

	///	set precision of numerical approximation of the tangent
		void set_tangent_precision(const number tanAccur)
		{m_tangentAccur = tanAccur;
		std::stringstream ss; ss << m_materialConfiguration << "accuracy of the tangent approximation = " << tanAccur << "\n";
		m_materialConfiguration = ss.str();}


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
			elasticityTensor(const size_t ip, MathMatrix<dim, dim>& GradU);

		virtual bool needs_to_add_jac_m(){return false;}


		virtual void init_internal_vars(TBaseElem* elem, const size_t numIP);
		virtual void internal_vars(TBaseElem* elem);
		virtual void update_internal_vars(const size_t ip, const MathMatrix<dim, dim>& GradU);

	///	use this method to make sure that all required attachments are attached
	/**	This method won't be necessary if we attach m_aElemData during initialization.*/
		virtual void attach_internal_vars(typename TDomain::grid_type& grid)
		{
	//todo:	Move this to an initialization routine. Best would probably to make
	//		set_domain(...) of the base class virtual and to attach m_aElemData there.

			grid.template attach_to<TBaseElem>(m_aElemData);
			m_aaElemData.access(grid, m_aElemData);
		}

		virtual void clear_attachments(typename TDomain::grid_type& grid)
		{
			if(grid.template has_attachment<TBaseElem>(m_aElemData)){
				grid.template detach_from<TBaseElem>(m_aElemData);
				m_aaElemData.invalidate();
			}
		}

	private:
		void strainTensor(MathMatrix<dim,dim>& strainTens, const MathMatrix<dim, dim>& GradU);

		void Flowrule(MathMatrix<dim, dim>& strain_p_new,
				MathMatrix<dim, dim>& strain,
				number& gamma,
				MathMatrix<dim, dim>& strial,
				MathMatrix<dim, dim>& normal,
				const MathMatrix<dim, dim>& GradU,
				const MathMatrix<dim, dim>& strain_p_old_t,
				const number alpha);

		void ConstLaw(MathMatrix<dim, dim>& stressTens, const MathMatrix<dim, dim>& strainTens,
				const MathMatrix<dim, dim>& strial, const number& gamma,
				const MathMatrix<dim, dim>& normal);

		void StressTensor(MathMatrix<dim,dim>& stressTens, const MathMatrix<dim, dim>& GradU,
						const MathMatrix<dim, dim>& strain_p_old_t, const number alpha);

		number Hardening(const number alpha);

		number Hardening_d(const number alpha);

		number PerfectPlasticity(const number flowcondtrial);

		number LinearHardening(const number flowcondtrial);

		number ExponentialHardening(const number strialnorm, const number alpha);

		void Update_internal_vars(MathMatrix<dim, dim>& strain_p_new,
				number& alpha,
				const MathMatrix<dim, dim>& GradU,
				const MathMatrix<dim, dim>& strain_p_old_t);

		//TODO: move this to common/math
	///	get trace and deviatoric part of a matrix
		number MatDeviatorTrace(const MathMatrix<dim, dim>& mat, MathMatrix<dim, dim>& dev);

	private:
	/// attached ElemData
		struct InternalVars {
			//	plastic part of strain tensor wrt to last time step
			MathMatrix<dim, dim> strain_p_old_t;

			//	hardening variable
			number alpha;
		};

	//	std-vector of InternalVars
		struct ElemData{
			std::vector<InternalVars> internalVars;
		};

		ElemData* m_pElemData;

	//	attachment type: attachment of ElemData
		typedef Attachment<ElemData> AElemData;
		AElemData m_aElemData;
	//	the attachment accessor
		typedef Grid::AttachmentAccessor<TBaseElem, AElemData>	ElemDataAccessor;
		ElemDataAccessor m_aaElemData;

	/// hardening behavior
		int m_hardening;
		size_t m_MaxHardIter;
		number m_HardAccuracy;

		struct MaterialConstants{
			number mu;		//	shear modulus
			number kappa;	//	bulk modulus
			number K_0; 	//	initial flow-stress:
			number K_inf; 	//	residual flow-stress (saturation stress)
			number Hard;	//	linear hardening modulus
			number omega;	//	hardening exponent
		}matConsts;

	/// tangent accuracy
		number m_tangentAccur;

	///	flags indicating if hardening variables are set
		bool m_bHardModulus;
		bool m_bHardExp;

	public:
		using base_type::m_materialConfiguration;
};

}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#include "prandtl_reuss_impl.h"

#endif /* PRANDTL_REUSS_H_ */
