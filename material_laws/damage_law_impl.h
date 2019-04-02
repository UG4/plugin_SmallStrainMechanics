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

#ifndef DAMAGE_LAW_IMPL_H_
#define DAMAGE_LAW_IMPL_H_

#include "damage_law.h"

#include "common/util/string_table_stream.h"

namespace ug{
namespace SmallStrainMechanics{

template <typename TDomain>
DamageLaw<TDomain>::
DamageLaw(	SmartPtr<GridFunction<TDomain,CPUAlgebra> > spPsi0,
			SmartPtr<GridFunction<TDomain,CPUAlgebra> > f)
	: IMaterialLaw<TDomain>(), 
	m_spPsi0(spPsi0), m_spF(f),
	m_pF_elem(NULL), m_pPsi0_elem(NULL)
{

}

template <typename TDomain>
void
DamageLaw<TDomain>::
init()
{
	//	check, if ElasticityTensor is set
	if (m_spElastTensorFunct.invalid())
		UG_THROW("No elasticity tensor set in DamageLaw::init()!");

	base_type::m_bInit = true;
}



template <typename TDomain>
void 
DamageLaw<TDomain>::
attach_internal_vars(typename TDomain::grid_type& grid)
{

}

template <typename TDomain>
void 
DamageLaw<TDomain>::
clear_attachments(typename TDomain::grid_type& grid)
{

}

template <typename TDomain>
void
DamageLaw<TDomain>::
init_internal_vars(TBaseElem* elem, const size_t numIP)
{
	if (!base_type::m_bInit){
		base_type::m_bInit = true;

		m_spF->set(1.0);
		m_spPsi0->set(0.0);
	}
}

template <typename TDomain>
inline
void
DamageLaw<TDomain>::
internal_vars(TBaseElem* elem)
{

	std::vector<DoFIndex> ind;
	const size_t fct = 0;

	if(m_spF->inner_dof_indices(elem, fct, ind) != 1)
		UG_THROW("Wrong number dofs");

	m_pF_elem = & DoFRef(*m_spF, ind[0]);
	m_pPsi0_elem = &DoFRef(*m_spPsi0, ind[0]);
}



template <typename TDomain>
void
DamageLaw<TDomain>::
stressTensor(MathMatrix<dim,dim>& stressTens, const size_t ip,
		const MathMatrix<dim, dim>& GradU)
{
	//	get linearized strain tensor (eps) at ip
	MathMatrix<dim, dim> strainTens;
	strainTensor(strainTens, GradU);

	//	TODO: replace this with general implementation of TensContractMat
	for(size_t i = 0; i < (size_t) dim; ++i)
		for(size_t j = 0; j < (size_t) dim; ++j)
		{
			stressTens[i][j] = 0.0;

			for(size_t k = 0; k < (size_t) dim; ++k)
				for(size_t l = 0; l < (size_t) dim; ++l)
					stressTens[i][j] += (*m_pF_elem) 
										* (*m_spElastTensorFunct)[i][j][k][l]
					                     * strainTens[k][l];

		}


}


template <typename TDomain>
SmartPtr<MathTensor4<TDomain::dim,TDomain::dim,TDomain::dim,TDomain::dim> >
DamageLaw<TDomain>::
elasticityTensor(const size_t ip, const MathMatrix<dim, dim>& GradU_const)
{
	return m_spElastTensorFunct;
}

template <typename TDomain>
inline
void
DamageLaw<TDomain>::
strainTensor(MathMatrix<dim,dim>& strainTens, const MathMatrix<dim, dim>& GradU)
{
	for(size_t i = 0; i < (size_t) dim; ++i)
		for(size_t j = 0; j < (size_t) dim; ++j)
			strainTens[i][j] = 0.5 * (GradU[i][j] + GradU[j][i]);
}

template <typename TDomain>
void
DamageLaw<TDomain>::
set_elasticity_tensor_orthotropic_plain_stress_E_G_nu(
const number E1, const number E2, //const number E3, 
const number G12, //const number G13, const number G23, 
const number v12 //const number v13, const number v23
)
{
	if(dim != 2) UG_THROW("Orthotrope Tensor PLAIN STRESS only for 2 dimensions" );

	MathTensor4<dim,dim,dim,dim> elastTensorFunct;

	//  setze alle wWerte auf 0
	for (size_t i = 0; i < (size_t) dim; ++i)
		for (size_t j = 0; j < (size_t) dim; ++j)
			for (size_t k = 0; k < (size_t) dim; ++k)
				for (size_t l = 0; l < (size_t) dim; ++l)
					elastTensorFunct[i][j][k][l] = 0.0;


	const number v21 = (E2/E1) * v12;
	//const number v31 = (E3/E1) * v13;
	//const number v32 = (E3/E2) * v23;

	const number D = 1 - v12*v21;


	// Tensor mit Werte fuellen
	//                 i  j  k  l
	elastTensorFunct[0][0][0][0] = E1/D; // C11
	elastTensorFunct[1][1][1][1] = E2/D; // C22

	elastTensorFunct[0][0][1][1] = 
	elastTensorFunct[1][1][0][0] = E2*v12/D; // = C12

	elastTensorFunct[0][1][0][1] = 
	elastTensorFunct[1][0][0][1] = 
	elastTensorFunct[0][1][1][0] = 
	elastTensorFunct[1][0][1][0] = 2*G12; // C66

	//	remembering the elasticity tensor
	SmartPtr<MathTensor4<dim,dim,dim,dim> > spElastTens(new MathTensor4<dim,dim,dim,dim>(elastTensorFunct));
	m_spElastTensorFunct = spElastTens;

	StringTableStream sts;
	std::stringstream ss;

	ss << "Hooke constants Orthotrope (PLAIN STRESS):\n";
	sts.clear();
	sts << " E1 = " << E1 << "E2 = " << E2 << "\n";
	sts << " v12 = " << v12 << "v21 = " << v21 << "\n";
	sts << "  G12 = " << G12 << "\n";
	ss << sts;
	m_materialConfiguration = ss.str();
	UG_LOG("\nset_elasticity_tensor_orthotropic " << ConfigShift(m_materialConfiguration) << "\n");
}


template <typename TDomain>
void
DamageLaw<TDomain>::
set_elasticity_tensor_orthotropic_plain_strain_E_G_nu(
const number E1, const number E2, const number E3, 
const number G12, const number G13, const number G23, 
const number v12, const number v13, const number v23
)
{
	if(dim != 2) UG_THROW("Orthotrope Tensor PLAIN STRAIN only for 2 dimensions" );

	MathTensor4<dim,dim,dim,dim> elastTensorFunct;

	//  setze alle wWerte auf 0
	for (size_t i = 0; i < (size_t) dim; ++i)
		for (size_t j = 0; j < (size_t) dim; ++j)
			for (size_t k = 0; k < (size_t) dim; ++k)
				for (size_t l = 0; l < (size_t) dim; ++l)
					elastTensorFunct[i][j][k][l] = 0.0;


	const number v21 = (E2/E1) * v12;
	const number v31 = (E3/E1) * v13;
	const number v32 = (E3/E2) * v23;

	const number D = 1 - v12*v21 - v13*v31 - v23*v32 - 2*v12*v23*v31;


	// Tensor mit Werte fuellen
	//                 i  j  k  l
	elastTensorFunct[0][0][0][0] = E1*(1-v23*v32)/D; // C11
	elastTensorFunct[1][1][1][1] = E2*(1-v13*v31)/D; // C22

	elastTensorFunct[0][0][1][1] = 
	elastTensorFunct[1][1][0][0] = E2*(v13*v32+v12)/D; // = C12

	elastTensorFunct[0][1][0][1] = 
	elastTensorFunct[1][0][0][1] = 
	elastTensorFunct[0][1][1][0] = 
	elastTensorFunct[1][0][1][0] = 2*G12; // C66

	//	remembering the elasticity tensor
	SmartPtr<MathTensor4<dim,dim,dim,dim> > spElastTens(new MathTensor4<dim,dim,dim,dim>(elastTensorFunct));
	m_spElastTensorFunct = spElastTens;

	StringTableStream sts;
	std::stringstream ss;

	ss << "Hooke constants Orthotrope (PLAIN STRAIN):\n";
	sts.clear();
	sts << " E1 = " << E1 << "E2 = " << E2 << "E3 = " << E3 << "\n";
	sts << " v12 = " << v12 << "v23 = " << v23 << "v13 = " << v13 << "\n";
	sts << "  G12 = " << G12 << "G23 = " << G23 << "G13 = " << G13 << "\n";
	ss << sts;
	m_materialConfiguration = ss.str();
	UG_LOG("\nset_elasticity_tensor_orthotropic " << ConfigShift(m_materialConfiguration) << "\n");
}



template <typename TDomain>
void
DamageLaw<TDomain>::
set_elasticity_tensor_orthotropic_E_G_nu(
const number E1, const number E2, const number E3, 
const number G12, const number G13, const number G23, 
const number v12, const number v13, const number v23
)
{
	if(dim != 3) UG_THROW("Orthotrope Tensor only for 3 dimensions" );

	MathTensor4<dim,dim,dim,dim> elastTensorFunct;

	//  setze alle wWerte auf 0
	for (size_t i = 0; i < (size_t) dim; ++i)
		for (size_t j = 0; j < (size_t) dim; ++j)
			for (size_t k = 0; k < (size_t) dim; ++k)
				for (size_t l = 0; l < (size_t) dim; ++l)
					elastTensorFunct[i][j][k][l] = 0.0;


	const number v21 = (E2/E1) * v12;
	const number v31 = (E3/E1) * v13;
	const number v32 = (E3/E2) * v23;

	const number D = 1 - v12*v21 - v13*v31 - v23*v32 - 2*v12*v23*v31;


	// Tensor mit Werte fuellen
	//                 i  j  k  l
	elastTensorFunct[0][0][0][0] = E1*(1-v23*v32)/D; // C11
	elastTensorFunct[1][1][1][1] = E2*(1-v13*v31)/D; // C22
	elastTensorFunct[2][2][2][2] = E3*(1-v12*v21)/D; // C33

	elastTensorFunct[0][0][1][1] = 
	elastTensorFunct[1][1][0][0] = E2*(v13*v32+v12)/D; // = C12

	elastTensorFunct[0][0][2][2] = 
	elastTensorFunct[2][2][0][0] = E3*(v12*v23+v13)/D; // = C13

	elastTensorFunct[1][1][2][2] = 
	elastTensorFunct[2][2][1][1] = E3*(v21*v13+v23)/D; // = C23


	elastTensorFunct[1][2][1][2] = 
	elastTensorFunct[1][2][2][1] = 
	elastTensorFunct[2][1][1][2] = 
	elastTensorFunct[2][1][2][1] = 2*G23; // C44

	elastTensorFunct[2][0][2][0] = 
	elastTensorFunct[0][2][2][0] = 
	elastTensorFunct[2][0][0][2] = 
	elastTensorFunct[0][2][0][2] = 2*G13; // C55

	elastTensorFunct[0][1][0][1] = 
	elastTensorFunct[1][0][0][1] = 
	elastTensorFunct[0][1][1][0] = 
	elastTensorFunct[1][0][1][0] = 2*G12; // C66

	//	remembering the elasticity tensor
	SmartPtr<MathTensor4<dim,dim,dim,dim> > spElastTens(new MathTensor4<dim,dim,dim,dim>(elastTensorFunct));
	m_spElastTensorFunct = spElastTens;

	StringTableStream sts;
	std::stringstream ss;

	ss << "Hooke constants Orthotrope:\n";
	sts.clear();
	sts << " E1 = " << E1 << "E2 = " << E2 << "E3 = " << E3 << "\n";
	sts << " v12 = " << v12 << "v23 = " << v23 << "v13 = " << v13 << "\n";
	sts << "  G12 = " << G12 << "G23 = " << G23 << "G13 = " << G13 << "\n";
	ss << sts;
	m_materialConfiguration = ss.str();
	UG_LOG("\nset_elasticity_tensor_orthotropic " << ConfigShift(m_materialConfiguration) << "\n");
}


template <typename TDomain>
void
DamageLaw<TDomain>::
set_elasticity_tensor_orthotropic(
const number C11, const number C12, const number C13,
		const number C22, const number C23,
					const number C33,
							const number C44,
									const number C55,
											const number C66 )
{
	if(dim != 3) UG_THROW("Orthotrope Tensor only for 3 dimensions" );

	MathTensor4<dim,dim,dim,dim> elastTensorFunct;

	//  setze alle wWerte auf 0
	for (size_t i = 0; i < (size_t) dim; ++i)
		for (size_t j = 0; j < (size_t) dim; ++j)
			for (size_t k = 0; k < (size_t) dim; ++k)
				for (size_t l = 0; l < (size_t) dim; ++l)
					elastTensorFunct[i][j][k][l] = 0.0;


	// Tensor mit Werte fuellen
	//                 i  j  k  l
	elastTensorFunct[0][0][0][0] = C11;

	elastTensorFunct[0][0][1][1] = C12;
	elastTensorFunct[1][1][0][0] = C12; // = C21

	elastTensorFunct[0][0][2][2] = C13;
	elastTensorFunct[2][2][0][0] = C13; // = C31

	elastTensorFunct[1][1][1][1] = C22;

	elastTensorFunct[1][1][2][2] = C23;
	elastTensorFunct[2][2][1][1] = C23; // = C32

	elastTensorFunct[2][2][2][2] = C33;

	elastTensorFunct[1][2][1][2] = C44;
	elastTensorFunct[1][2][2][1] = C44;
	elastTensorFunct[2][1][1][2] = C44;
	elastTensorFunct[2][1][2][1] = C44;

	elastTensorFunct[2][0][2][0] = C55;
	elastTensorFunct[0][2][2][0] = C55;
	elastTensorFunct[2][0][0][2] = C55;
	elastTensorFunct[0][2][0][2] = C55;

	elastTensorFunct[0][1][0][1] = C66;
	elastTensorFunct[1][0][0][1] = C66;
	elastTensorFunct[0][1][1][0] = C66;
	elastTensorFunct[1][0][1][0] = C66;

	//	remembering the elasticity tensor
	SmartPtr<MathTensor4<dim,dim,dim,dim> > spElastTens(new MathTensor4<dim,dim,dim,dim>(elastTensorFunct));
	m_spElastTensorFunct = spElastTens;

	DenseMatrix<FixedArray2<double, 3, 3> > mat;
	mat(0,0) = C11; mat(0,1) = C12; mat(0,2) = C13;
	mat(1,0) = C12; mat(1,1) = C22; mat(1,2) = C23;
	mat(2,0) = C13; mat(2,1) = C23; mat(2,2) = C33;

	Invert(mat);
	//UG_LOG("S = " mat << "\n");
	number E1 = 1/mat(0,0), E2 = 1/mat(1,1), E3 = 1/mat(2,2);
	number v12 = -mat(1, 0)/E1;
	//number v21 = -mat(0, 1)/E2;
	number v23 = -mat(2, 1)/E2;
	//number v32 = -mat(1, 2)/E3;
	number v13 = -mat(2, 0)/E1;
	//number v31 = -mat(0, 2)/E3;
	number G23 = C44;
	number G13 = C55;
	number G12 = C66;

	StringTableStream sts;
	std::stringstream ss;
	ss << "Orthrope Elasticity Tensor: \n";
	sts.clear();
	sts << C11 << C12 << C13 << "\n";
	sts << C12 << C22 << C23 << "\n";
	sts << C13 << C23 << C33 << "\n";
	sts << sts.empty_col(3) << C44 << "\n";
	sts << sts.empty_col(4)  << C55 << "\n";
	sts << sts.empty_col(5)  << C66 << "\n";
	ss << sts;

	ss << "Hooke constants Orthotrope:\n";
	sts.clear();
	sts << " E1 = " << E1 << "E2 = " << E2 << "E3 = " << E3 << "\n";
	sts << " v12 = " << v12 << "v23 = " << v23 << "v13 = " << v13 << "\n";
	sts << "  G12 = " << G12 << "G23 = " << G23 << "G13 = " << G13 << "\n";
	ss << sts;
	m_materialConfiguration = ss.str();
	UG_LOG("\nset_elasticity_tensor_orthotropic " << ConfigShift(m_materialConfiguration) << "\n");
}


template <typename TDomain>
void
DamageLaw<TDomain>::
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
DamageLaw<TDomain>::
set_hooke_elasticity_tensor_plain_strain_E_nu(const number E, const number nu)
{
	if(dim != 2) UG_THROW("Plain Strain Tensor only for 2 dimensions" );

	set_hooke_elasticity_tensor_E_nu(E, nu);
}


template <typename TDomain>
void
DamageLaw<TDomain>::
set_hooke_elasticity_tensor_E_nu(const number E, const number nu)
{
	number lambda = (E*nu) / ((1+nu)*(1-2*nu));
	number mu = E/(2*(1+nu));
	set_hooke_elasticity_tensor(lambda, mu);
}

template <typename TDomain>
void
DamageLaw<TDomain>::
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

}
}

#endif /* DAMAGE_LAW_IMPL_H_ */
