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

#ifndef SKIN_LAW_IMPL_H_
#define SKIN_LAW_IMPL_H_

#include "common/util/string_table_stream.h"

#include "skin_law.h"

#define PROFILE_HOOKE
#ifdef PROFILE_HOOKE
	#define HOOKE_PROFILE_FUNC()		PROFILE_FUNC_GROUP("Hooke")
	#define HOOKE_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "Hooke")
	#define HOOKE_PROFILE_END()		PROFILE_END()
#else
	#define HOOKE_PROFILE_FUNC()
	#define HOOKE_PROFILE_BEGIN(name)
	#define HOOKE_PROFILE_END()
#endif

namespace ug{
namespace SmallStrainMechanics{

template <typename TDomain>
void
SkinMaterialLaw<TDomain>::
init()
{
	// dummy values
	set_hooke_elasticity_tensor(-10000, -100000);

	//	check, if ElasticityTensor is set
	if (m_spElastTensorFunct.invalid())
		UG_THROW("No elasticity tensor set in SkinMaterialLaw::init()!");

	base_type::m_bInit = true;
}

template <typename TDomain>
inline
void
SkinMaterialLaw<TDomain>::
strainTensor(MathMatrix<dim,dim>& strainTens, const MathMatrix<dim, dim>& GradU)
{
	for(size_t i = 0; i < (size_t) dim; ++i)
		for(size_t j = 0; j < (size_t) dim; ++j)
			strainTens[i][j] = 0.5 * (GradU[i][j] + GradU[j][i]);
}


template <typename TDomain>
void
SkinMaterialLaw<TDomain>::
stressTensor(MathMatrix<dim,dim>& stressTens, const size_t ip,
		const MathVector<dim>& x, const MathMatrix<dim, dim>& GradU)
{
	HOOKE_PROFILE_BEGIN(SkinMaterialLaw_stressTensor);

	const number y = x[1];
	number _y = 1.0;

	if      (0.9375 >= y && y > 0.8125)	_y = 0.875;
	else if (0.8125 >= y && y > 0.6875) _y = 0.75;
	else if (0.6875 >= y && y > 0.5625) _y = 0.625;
	else if (0.5625 >= y && y > 0.4375) _y = 0.5;
	else if (0.4375 >= y && y > 0.3125) _y = 0.375;
	else if (0.3125 >= y && y > 0.1875) _y = 0.25;
	else if (0.1875 >= y && y > 0.0625) _y = 0.125;

	if(_y == 1.0) UG_THROW("Position not found: " << y);

	const number lambda_top = 1.0, lambda_bottom = 2.0;
	const number mu_top = 0.1, mu_bottom = 0.2;

	const number mu = mu_bottom * (1.0 - _y) + mu_top * (_y);
	const number lambda = lambda_bottom * (1.0 - _y) + lambda_top * (_y);


/*
	const number lambda_top = 1.0, lambda_bottom = 2.0;
	const number mu_top = 0.1, mu_bottom = 0.2;

	const number mu = mu_bottom * (1.0 - x[1]) + mu_top * (x[1]);
	const number lambda = lambda_bottom * (1.0 - x[1]) + lambda_top * (x[1]);
*/

	//	get linearized strain tensor (eps) at ip
	MathMatrix<dim, dim> strainTens;
	strainTensor(strainTens, GradU);

	MatScale(stressTens, 2*mu, strainTens);

	const number trace = Trace(strainTens);
	for (int d = 0; d < dim; ++d)
		stressTens[d][d] += lambda * trace; 

}

template <typename TDomain>
SmartPtr<MathTensor4<TDomain::dim,TDomain::dim,TDomain::dim,TDomain::dim> >
SkinMaterialLaw<TDomain>::
elasticityTensor(const size_t ip, const MathVector<dim>& x, const MathMatrix<dim, dim>& GradU)
{
	const number y = x[1];
	number _y = 1.0;

	if      (0.9375 >= y && y > 0.8125)	_y = 0.875;
	else if (0.8125 >= y && y > 0.6875) _y = 0.75;
	else if (0.6875 >= y && y > 0.5625) _y = 0.625;
	else if (0.5625 >= y && y > 0.4375) _y = 0.5;
	else if (0.4375 >= y && y > 0.3125) _y = 0.375;
	else if (0.3125 >= y && y > 0.1875) _y = 0.25;
	else if (0.1875 >= y && y > 0.0625) _y = 0.125;

	if(_y == 1.0) UG_THROW("Position not found: " << y);

	const number lambda_top = 1.0, lambda_bottom = 2.0;
	const number mu_top = 0.1, mu_bottom = 0.2;

	const number mu = mu_bottom * (1.0 - _y) + mu_top * (_y);
	const number lambda = lambda_bottom * (1.0 - _y) + lambda_top * (_y);

	set_hooke_elasticity_tensor(lambda, mu);

	return m_spElastTensorFunct;
}


template <typename TDomain>
void
SkinMaterialLaw<TDomain>::
set_hooke_elasticity_tensor_plain_stress_E_nu(const number E, const number nu)
{
	if(dim != 2) UG_THROW("Plain Stress Tensor only for 2 dimensions" );

	MathTensor4<dim,dim,dim,dim> elastTensorFunct;

	//  setze alle wWerte auf 0
	for (size_t i = 0; i < (size_t) dim; ++i)
		for (size_t j = 0; j < (size_t) dim; ++j)
			for (size_t k = 0; k < (size_t) dim; ++k)
				for (size_t l = 0; l < (size_t) dim; ++l)
					elastTensorFunct[i][j][k][l] = 0.0;

	// Tensor mit Werte fuellen
	//                 i  j  k  l
	elastTensorFunct[0][0][0][0] = E/(1-nu*nu); // C11
	elastTensorFunct[1][1][1][1] = E/(1-nu*nu); // C22

	elastTensorFunct[0][0][1][1] = 
	elastTensorFunct[1][1][0][0] = (E*nu)/(1-nu*nu); // = C12 = C21

	elastTensorFunct[0][1][0][1] = 
	elastTensorFunct[1][0][0][1] = 
	elastTensorFunct[0][1][1][0] = 
	elastTensorFunct[1][0][1][0] = E*(1-nu)/(1-nu*nu); // C33

	//	remembering the elasticity tensor
	SmartPtr<MathTensor4<dim,dim,dim,dim> > spElastTens(new MathTensor4<dim,dim,dim,dim>(elastTensorFunct));
	m_spElastTensorFunct = spElastTens;

	std::stringstream ss;
	ss << "Hooke Elasticity Isotropic (PLAIN STRESS): \n";
	ss << "  young modulus (Elastizitaetsmodul): " << E << "\n";
	ss << "  poisson ratio (Querkontraktionszahl) v: " << nu << "\n";
	ss << "  Elasticity Tensor = " << elastTensorFunct << "\n";
	m_materialConfiguration = ss.str();
}

template <typename TDomain>
void
SkinMaterialLaw<TDomain>::
set_hooke_elasticity_tensor_plain_strain_E_nu(const number E, const number nu)
{
	if(dim != 2) UG_THROW("Plain Strain Tensor only for 2 dimensions" );

	set_hooke_elasticity_tensor_E_nu(E, nu);
}


template <typename TDomain>
void
SkinMaterialLaw<TDomain>::
set_hooke_elasticity_tensor_E_nu(const number E, const number nu)
{
	number lambda = (E*nu) / ((1+nu)*(1-2*nu));
	number mu = E/(2*(1+nu));
	set_hooke_elasticity_tensor(lambda, mu);
}

template <typename TDomain>
void
SkinMaterialLaw<TDomain>::
set_hooke_elasticity_tensor(const number lambda, const number mu)
{
	//	sets the 'Hooke'-elasticity tensor for isotropic materials
	//	in 2D this tensor formulation corresponds to the
	//	plane strain assumption for Hookes`s law

	MathTensor4<dim,dim,dim,dim> elastTensorFunct;

	//  filling the constant elasticity tensor
	for (size_t i = 0; i < (size_t) dim; ++i)
		for (size_t j = 0; j < (size_t) dim; ++j)
			for (size_t k = 0; k < (size_t) dim; ++k)
				for (size_t l = 0; l < (size_t) dim; ++l)
				{
					elastTensorFunct[i][j][k][l] = 0.0;

					if ((i == j) && (k == l))
						elastTensorFunct[i][j][k][l] += lambda;

					if ((i == k) && (j == l))
						elastTensorFunct[i][j][k][l] +=  mu;

					if ((i == l) && (j == k))
						elastTensorFunct[i][j][k][l] +=  mu;

				} //end (l)

	//	remembering the elasticity tensor
	SmartPtr<MathTensor4<dim,dim,dim,dim> > spElastTens(new MathTensor4<dim,dim,dim,dim>(elastTensorFunct));
	m_spElastTensorFunct = spElastTens;

	number E = mu * (3.0 * lambda + 2.0 * mu)/(lambda + mu);
	number v = 0.5 * lambda/(lambda + mu);
	number kappa = lambda + 2.0/3.0 * mu;

	std::stringstream ss;
	ss << "Hooke Elasticity Isotropic: \n";
	ss << "  Lame`s first constant lambda: " << lambda << "\n";
	ss << "  Lame`s second constant mue (sometimes 'G', shear modulus) (Schubmodul): " << mu << "\n";
	ss << " This setting equals: \n";
	ss << "  young modulus (Elastizitaetsmodul): " << E << "\n";
	ss << "  poisson ratio (Querkontraktionszahl) v: " << v << "\n";
	ss << "  bulk modulus (Kompressionsmodul): " << kappa << "\n";
	ss << "  Elasticity Tensor = " << elastTensorFunct << "\n";
	m_materialConfiguration = ss.str();
}

}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#endif /* SKIN_LAW_IMPL_H_ */
