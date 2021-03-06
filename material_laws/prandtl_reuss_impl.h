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

#ifndef PRANDTL_REUSS_IMPL_H_
#define PRANDTL_REUSS_IMPL_H_

#include "prandtl_reuss.h"

#define PROFILE_PRANDTL_REUSS
#ifdef PROFILE_PRANDTL_REUSS
	#define PRANDTL_REUSS_PROFILE_FUNC()		PROFILE_FUNC_GROUP("Prandtl Reuss")
	#define PRANDTL_REUSS_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "Prandtl Reuss")
	#define PRANDTL_REUSS_PROFILE_END()		PROFILE_END()
#else
	#define PRANDTL_REUSS_PROFILE_FUNC()
	#define PRANDTL_REUSS_PROFILE_BEGIN(name)
	#define PRANDTL_REUSS_PROFILE_END()
#endif

namespace ug{
namespace SmallStrainMechanics{

template <typename TDomain>
PrandtlReuss<TDomain>::PrandtlReuss():
	IMaterialLaw<TDomain>(),
	m_MaxHardIter(100), m_HardAccuracy(0.0), m_tangentAccur(1e-08),
	m_bHardModulus(false), m_bHardExp(false)
{
	// set default material constants
	matConsts.mu = 0.0; matConsts.kappa = 0.0;

	matConsts.K_0 = 0.0; m_hardening = 0;
	matConsts.K_inf = 0.0; matConsts.Hard = 0.0; matConsts.omega = 0.0;

	m_max_k_tan = 0.0;
	m_min_k_tan = 100.0;
	m_plasticIPs = 0;

	std::stringstream ss;
	ss << "Prandtl Reuss Plasticity: \n";
	m_materialConfiguration = ss.str();

	base_type::m_bConstElastTens = false;
}

template <typename TDomain>
void
PrandtlReuss<TDomain>::
set_hardening_behavior(int hard)
{
	std::stringstream ss;

	switch (hard){
		//	perfect hardening
		case 0: m_hardening = 0;
				ss << m_materialConfiguration << "perfect hardening \n";
				m_materialConfiguration = ss.str();
				break;

		//	linear hardening
		case 1: m_hardening = 1;
				ss << m_materialConfiguration << "linear hardening \n";
				m_materialConfiguration = ss.str();
				break;

		//	exponential hardening
		case 2: m_hardening = 2;
				m_MaxHardIter = 100;
				m_HardAccuracy = 1e-10;
				ss << m_materialConfiguration << "exponential hardening \n"
				<< " max. hardening iterations = " << m_MaxHardIter << "\n"
				<< " hardening accuracy = " << m_HardAccuracy << "\n";
				m_materialConfiguration = ss.str();
				break;

		default: UG_THROW(hard << " is not a valid hardening behavior! "
				"Choose 0 (perfect), 1 (linear) or 2 (exponential) ! \n");
	}
}

template <typename TDomain>
void
PrandtlReuss<TDomain>::
init()
{
	//	check, if all parameter for a hardening behavior are set
	if (m_hardening == 1){ // linear hardening law
		if (!m_bHardModulus)
			UG_THROW("No hardening modulus set! This is necessary for a "
					"linear hardening law in PrandtlReuss::init() \n");
	}

	if (m_hardening == 2){ // exponential hardening law
		if (!m_bHardModulus)
			UG_THROW("No hardening modulus set! This is necessary for an "
					"exponential hardening law in PrandtlReuss::init() \n");
		if (!m_bHardExp)
			UG_THROW("No hardening exponent set! This is necessary for an "
					"exponential hardening law in PrandtlReuss::init() \n");
	}
}

template <typename TDomain>
void
PrandtlReuss<TDomain>::
StressTensor(MathMatrix<dim,dim>& stressTens, const MathMatrix<dim, dim>& GradU,
		const MathMatrix<dim, dim>& strain_p_old_t, const number alpha)
{
	MathMatrix<dim, dim> strial, normal, strain, strain_p_new;
	number gamma;

	Flowrule(strain_p_new, strain, gamma, strial, normal, GradU, strain_p_old_t, alpha);
	ConstLaw(stressTens, strain, strial, gamma, normal);
}

template <typename TDomain>
void
PrandtlReuss<TDomain>::
stressTensor(MathMatrix<dim,dim>& stressTens, const size_t ip,
		const MathMatrix<dim, dim>& GradU)
{
	PRANDTL_REUSS_PROFILE_BEGIN(PrandtlReuss_stressTensor);

	//	get internal variables
	MathMatrix<dim, dim>& strain_p_old_t = m_pElemData->internalVars[ip].strain_p_old_t;
	number alpha = m_pElemData->internalVars[ip].alpha;

	StressTensor(stressTens, GradU, strain_p_old_t, alpha);
}


template <typename TDomain>
SmartPtr<MathTensor4<TDomain::dim,TDomain::dim,TDomain::dim,TDomain::dim> >
PrandtlReuss<TDomain>::
elasticityTensor(const size_t ip, const MathMatrix<dim, dim>& GradU_const)
{
	PRANDTL_REUSS_PROFILE_BEGIN(PrandtlReuss_elasticityTensor);

	//	get internal variables
	MathMatrix<dim, dim>& strain_p_old_t = m_pElemData->internalVars[ip].strain_p_old_t;
	number alpha = m_pElemData->internalVars[ip].alpha;

	MathMatrix<dim,dim>& GradU = const_cast<MathMatrix<dim,dim>&>(GradU_const);

	//	computes the elasticity Tensor by means of numerical differentiation
	MathTensor4<dim,dim,dim,dim> elastTens;
	MathMatrix<dim, dim> stressT, stressTT;

	//	for formulation in reference config
	StressTensor(stressT, GradU, strain_p_old_t, alpha);

	for (size_t i = 0; i < (size_t) dim; ++i)
		for (size_t j = 0; j < (size_t) dim; ++j)
		{
			GradU[i][j] += m_tangentAccur;

			StressTensor(stressTT, GradU, strain_p_old_t, alpha);

			for (size_t k = 0; k < (size_t) dim; ++k)
				for (size_t l = 0; l < (size_t) dim; ++l)
					elastTens[i][j][k][l] = 1.0/ m_tangentAccur * (stressTT[k][l] - stressT[k][l]);

			GradU[i][j] -= m_tangentAccur;
		}

	//	compute maximal relative error of numerical differentiation
	//	(ref: Deuflhard Numerische Mathematik 1)
	MathMatrix<dim, dim> diffStress;
	MatSubtract(diffStress, stressTT, stressT);
	if (MatFrobeniusNorm(stressT) > 0.0)
	{
		const number k_tan = MatFrobeniusNorm(diffStress)/MatFrobeniusNorm(stressT);
		if (k_tan > m_max_k_tan)
			m_max_k_tan = k_tan;
		if (k_tan < m_min_k_tan)
			m_min_k_tan = k_tan;
	}

	//	TODO: change this smart pointer to a member variable
	//	do the same with m_GradU!
	SmartPtr<MathTensor4<dim,dim,dim,dim> > spElastTens(new MathTensor4<dim,dim,dim,dim>(elastTens));
	return spElastTens;
}

template <typename TDomain>
void
PrandtlReuss<TDomain>::
init_internal_vars(TBaseElem* elem, const size_t numIP)
{
	m_aaElemData[elem].internalVars.resize(numIP);

	// 	set plastic strain (eps_p) and hardening variable (alpha)
	//	to zero at every ip (in ElemData-struct)
	for (size_t ip = 0; ip < numIP; ++ip)
	{
		m_aaElemData[elem].internalVars[ip].strain_p_old_t = 0.0;
		m_aaElemData[elem].internalVars[ip].alpha = 0.0;
	}

	if (!base_type::m_bInit)
		base_type::m_bInit = true;
}

template <typename TDomain>
inline
void
PrandtlReuss<TDomain>::
internal_vars(TBaseElem* elem)
{
	m_pElemData = &m_aaElemData[elem];
}

template <typename TDomain>
void
PrandtlReuss<TDomain>::
update_internal_vars(const size_t ip, const MathMatrix<dim, dim>& GradU)
{
	MathMatrix<dim, dim>& strain_p_old_t = m_pElemData->internalVars[ip].strain_p_old_t;
	//	TODO use a reference for alpha here!
	number alpha = m_pElemData->internalVars[ip].alpha;

	//  update the internal variables: strain_p_old_t (plastic strain tensor)
	//	and alpha (hardening parameter) at ip
	Update_internal_vars(strain_p_old_t, alpha, GradU, strain_p_old_t);
	m_pElemData->internalVars[ip].alpha = alpha;
}

template <typename TDomain>
void
PrandtlReuss<TDomain>::
Update_internal_vars(MathMatrix<dim, dim>& strain_p_new,
		number& alpha,
		const MathMatrix<dim, dim>& GradU,
		const MathMatrix<dim, dim>& strain_p_old_t)
{
	//	compute trial strain tensor (eps)
	MathMatrix<dim, dim> strainTensTrial;
	for(size_t i = 0; i < (size_t) dim; ++i)
		for(size_t j = 0; j < (size_t) dim; ++j)
			strainTensTrial[i][j] = 0.5 * (GradU[i][j] + GradU[j][i]) - strain_p_old_t[i][j];

	MathMatrix<dim, dim> dev_strainTrial;
	MatDeviatorTrace(strainTensTrial, dev_strainTrial);

	MathMatrix<dim, dim> strial;
	for(size_t i = 0; i < (size_t) dim; ++i)
		for(size_t j = 0; j < (size_t) dim; ++j)
			strial[i][j] = 2.0 * matConsts.mu * dev_strainTrial[i][j];

	number strialnorm = MatFrobeniusNorm(strial);
	number flowcondtrial = strialnorm - sqrt(2.0 / 3.0) * Hardening(alpha);

	if (flowcondtrial <= 0){
		strain_p_new = strain_p_old_t; return;
	}

	m_plasticIPs += 1;

	number gamma = 0.0;

	UG_ASSERT(strialnorm > 0.0, "norm of strial needs to be > 0.0");

	//	computation of gamma (plastic corrector/multiplicator)
	switch (m_hardening)
	{
		case 0: gamma = PerfectPlasticity(flowcondtrial); break;
		case 1: gamma = LinearHardening(flowcondtrial); break;
		case 2: gamma = ExponentialHardening(strialnorm, alpha); break;
		default:
			UG_THROW(m_hardening << " in 'Update' is not a valid hardening behavior! \n");
	}

	MathMatrix<dim, dim> normaltrial;
	MatScale(normaltrial, 1.0 / strialnorm, strial);

	UG_ASSERT(gamma > 0.0, "gamma needs to be > 0.0");

	MathMatrix<dim, dim> d_strain_p;
	MatScale(d_strain_p, gamma, normaltrial);
	MatAdd(strain_p_new, strain_p_old_t, d_strain_p);

	alpha += sqrt(2.0 / 3.0) * gamma;
}

template <typename TDomain>
inline
void
PrandtlReuss<TDomain>::
strainTensor(MathMatrix<dim,dim>& strainTens, const MathMatrix<dim, dim>& GradU)
{
	for(size_t i = 0; i < (size_t) dim; ++i)
		for(size_t j = 0; j < (size_t) dim; ++j)
			strainTens[i][j] = 0.5 * (GradU[i][j] + GradU[j][i]);
}

template <typename TDomain>
void
PrandtlReuss<TDomain>::
Flowrule(MathMatrix<dim, dim>& strain_p_new, MathMatrix<dim, dim>& strain, number& gamma,
		MathMatrix<dim, dim>& strial, MathMatrix<dim, dim>& normal, const MathMatrix<dim, dim>& GradU,
		const MathMatrix<dim, dim>& strain_p_old_t, const number alpha)
{
	//////////////////////////
	//  TRIAL ELASTIC STEP
	//////////////////////////

	MathMatrix<dim, dim> strain_trial;

	//	compute linearized strain tensor (eps)
	for(size_t i = 0; i < (size_t) dim; ++i)
		for(size_t j = 0; j < (size_t) dim; ++j)
		{
			strain[i][j] = 0.5 * (GradU[i][j] + GradU[j][i]);
			//	get eps_trial := eps_{n+1} - eps^p_{n}
			strain_trial[i][j] = strain[i][j] - strain_p_old_t[i][j];
		}

	MathMatrix<dim, dim> dev_strain_trial;
	MatDeviatorTrace(strain_trial, dev_strain_trial);

	//	compute trial strain deviator
	MatScale(strial, 2.0 * matConsts.mu, dev_strain_trial);

	number strialnorm = MatFrobeniusNorm(strial);
	number flowcondtrial = strialnorm - sqrt(2.0 / 3.0) * Hardening(alpha);

	//////////////////////////////
	//  CHECKING YIELD-CONDITION:
	//////////////////////////////

	if (flowcondtrial <= 0){
		strain_p_new = strain_p_old_t; gamma = 0.0;
		return;
	}

	////////////////////////////////////
	//  RETURN-MAPPING (corrector-step)
	////////////////////////////////////

	UG_ASSERT(strialnorm > 0.0, "norm of strial needs to be > 0.0");

	//	computation of gamma (plastic corrector/multiplicator)
	// 	accordingly to Simo/Hughes 98 p.121/122
	switch (m_hardening)
	{
		case 0: gamma = PerfectPlasticity(flowcondtrial); break;
		case 1: gamma = LinearHardening(flowcondtrial); break;
		case 2: gamma = ExponentialHardening(strialnorm, alpha); break;
		default:
			UG_THROW(m_hardening << " in 'Flowrule' is not a valid hardening behavior! \n");
	}

	UG_ASSERT(gamma > 0.0, "gamma needs to be > 0.0");

	MatScale(normal, 1.0 / strialnorm, strial);

	MathMatrix<dim, dim> d_strain_p;
	MatScale(d_strain_p, gamma, normal);
	MatAdd(strain_p_new, strain_p_old_t, d_strain_p);
}

template <typename TDomain>
void
PrandtlReuss<TDomain>::
ConstLaw(MathMatrix<dim, dim>& stressTens, const MathMatrix<dim, dim>& strain, const MathMatrix<dim, dim>& strial,
		const number& gamma, const MathMatrix<dim, dim>& normal)
{
	//	constitutive law taken from Simo/Hughes 98 'Computational Inelasticity' p. 124
	number trStrain = Trace(strain);

	//	compute sigma = kappa * tr[eps] * id + strial - 2 * mu * gamma * normal
	for(size_t i = 0; i < (size_t) dim; ++i)
	{
		for(size_t j = 0; j < (size_t) dim; ++j)
			stressTens[i][j] = strial[i][j] - 2.0 * matConsts.mu * gamma * normal[i][j];

		stressTens[i][i] += matConsts.kappa * trStrain;
	}

}

template <typename TDomain>
number
PrandtlReuss<TDomain>::
plastic_multiplier(const size_t ip, const MathMatrix<dim, dim>& GradU)
{
	//	get internal variables
	MathMatrix<dim, dim>& strain_p_old_t = m_pElemData->internalVars[ip].strain_p_old_t;
	number alpha = m_pElemData->internalVars[ip].alpha;

	/////////////////////////
	//  TRIAL ELASTIC STEP
	//////////////////////////

	MathMatrix<dim, dim> strain_trial;

	//	compute linearized strain tensor (eps)
	for(size_t i = 0; i < (size_t) dim; ++i)
		for(size_t j = 0; j < (size_t) dim; ++j)
		{
			//	get eps_trial := eps_{n+1} - eps^p_{n}
			strain_trial[i][j] = 0.5 * (GradU[i][j] + GradU[j][i]) - strain_p_old_t[i][j];
		}

	MathMatrix<dim, dim> dev_strain_trial;
	MatDeviatorTrace(strain_trial, dev_strain_trial);

	//	compute trial strain deviator
	MathMatrix<dim, dim> strial;
	MatScale(strial, 2.0 * matConsts.mu, dev_strain_trial);

	number strialnorm = MatFrobeniusNorm(strial);
	number flowcondtrial = strialnorm - sqrt(2.0 / 3.0) * Hardening(alpha);

	//////////////////////////////
	//  CHECKING YIELD-CONDITION:
	//////////////////////////////

	if (flowcondtrial <= 0){
		return 0.0;
	}

	////////////////////////////////////
	//  RETURN-MAPPING (corrector-step)
	////////////////////////////////////

	UG_ASSERT(strialnorm > 0.0, "norm of strial needs to be > 0.0");

	number gamma = 0.0;

	//	computation of gamma (plastic corrector/multiplicator)
	// 	accordingly to Simo/Hughes 98 p.121/122
	switch (m_hardening)
	{
		case 0: gamma = PerfectPlasticity(flowcondtrial); break;
		case 1: gamma = LinearHardening(flowcondtrial); break;
		case 2: gamma = ExponentialHardening(strialnorm, alpha); break;
		default:
			UG_THROW(m_hardening << " in 'Flowrule' is not a valid hardening behavior! \n");
	}

	UG_ASSERT(gamma > 0.0, "gamma needs to be > 0.0");

	return gamma;
}

template <typename TDomain>
inline
number
PrandtlReuss<TDomain>::
Hardening(const number alpha)
{
	return (matConsts.K_0 + matConsts.Hard * alpha +
			(matConsts.K_inf - matConsts.K_0) * (1.0 - exp(-matConsts.omega * alpha)));
}

//	"Hardening_d": derivative of nonlinear hardening law "Hardening"
template <typename TDomain>
inline
number
PrandtlReuss<TDomain>::
Hardening_d(const number alpha)
{
	return (matConsts.Hard + (matConsts.K_inf - matConsts.K_0) * matConsts.omega
			* exp(-matConsts.omega * alpha));
}

template <typename TDomain>
inline
number
PrandtlReuss<TDomain>::
PerfectPlasticity(const number flowcondtrial)
{
	return flowcondtrial / (2.0 * matConsts.mu);
}

template <typename TDomain>
inline
number
PrandtlReuss<TDomain>::
LinearHardening(const number flowcondtrial)
{
	return flowcondtrial / (2.0 * (matConsts.mu + matConsts.Hard/3.0));
}

template <typename TDomain>
number
PrandtlReuss<TDomain>::
ExponentialHardening(const number strialnorm, const number alpha)
{
	number gamma = 0.0;

	for (size_t i = 0; i < m_MaxHardIter; ++i) {
		//  f_cap has to be equal to 0, due to the consistency condition
		number f_cap = strialnorm - 2.0 * gamma * matConsts.mu - sqrt(2.0 / 3.0) * Hardening(
				alpha + sqrt(2.0 / 3.0) * gamma);

		//abs(f_cap) < HardAccuracy * strialnorm vs. - f_cap > - HardAccuracy * strialnorm
		if (f_cap < m_HardAccuracy * strialnorm){ break;}

		gamma += f_cap / (2.0 * matConsts.mu * (1.0 + Hardening_d(
				alpha + sqrt(2.0 / 3.0) * gamma) / (3.0 * matConsts.mu)));
	}

	return gamma;
}

template <typename TDomain>
void
PrandtlReuss<TDomain>::
write_data_to_console(const number t)
{
	UG_LOG("maximal k_tan:" << m_max_k_tan << "\n");
	UG_LOG("minimal k_tan:" << m_min_k_tan << "\n");

	//  print: at how many gauss points we are in the plastic zone,...
	UG_LOG("# of plastic IPs in this timestep:" << m_plasticIPs << "\n");
	//	reset data
	m_plasticIPs = 0;
}

}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#endif /* PRANDTL_REUSS_IMPL_H_ */
