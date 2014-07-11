/*
 * hooke_impl.h
 *
 *  Created on: 03.02.2014
 *      Author: raphaelprohl
 */

#ifndef HOOKE_IMPL_H_
#define HOOKE_IMPL_H_

#include "common/util/string_table_stream.h"

#include "hooke.h"

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
HookeLaw<TDomain>::
init()
{
	//	check, if ElasticityTensor is set
	if (m_spElastTensorFunct.invalid())
		UG_THROW("No elasticity tensor set in HookeLaw::init()!");

	base_type::m_bInit = true;
}

template <typename TDomain>
void
HookeLaw<TDomain>::
stressTensor(MathMatrix<dim,dim>& stressTens, const size_t ip,
		const MathMatrix<dim, dim>& GradU)
{
	HOOKE_PROFILE_BEGIN(HookeLaw_stressTensor);

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
					stressTens[i][j] += (*m_spElastTensorFunct)[i][j][k][l]
					                     * strainTens[k][l];

		}

}

template <typename TDomain>
inline
SmartPtr<MathTensor4<TDomain::dim,TDomain::dim,TDomain::dim,TDomain::dim> >
HookeLaw<TDomain>::
elasticityTensor()
{
	HOOKE_PROFILE_BEGIN(HookeLaw_elasticityTensor);
	return m_spElastTensorFunct;
}

template <typename TDomain>
inline
void
HookeLaw<TDomain>::
strainTensor(MathMatrix<dim,dim>& strainTens, const MathMatrix<dim, dim>& GradU)
{
	for(size_t i = 0; i < (size_t) dim; ++i)
		for(size_t j = 0; j < (size_t) dim; ++j)
			strainTens[i][j] = 0.5 * (GradU[i][j] + GradU[j][i]);
}

template <typename TDomain>
void
HookeLaw<TDomain>::
set_elasticity_tensor_orthotropic(
const number C11, const number C12, const number C13,
		const number C22, const number C23,
					const number C33,
							const number C44,
									const number C55,
											const number C66 )
{
	UG_ASSERT( dim==3, "Orthotrope Tensor only for 3 dimensions" );

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

	ss << "Hooke constants:\n";
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
HookeLaw<TDomain>::
set_hooke_elasticity_tensor_E_nu(const number E, const number nu)
{
	number lambda = (E*nu) / ((1+nu)*(1-2*nu));
	number mu = E/(2*(1+nu));
	set_hooke_elasticity_tensor(lambda, mu);
}

template <typename TDomain>
void
HookeLaw<TDomain>::
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
	ss << "Hooke Elasticity Tensor: \n";
	ss << "  Lame`s first constant lambda: " << lambda << "\n";
	ss << "  Lame`s second constant mue (sometimes 'G', shear modulus) (Schubmodul): " << mu << "\n";
	ss << " This setting equals: \n";
	ss << "  young modulus (Elastizitaetsmodul): " << E << "\n";
	ss << "  poisson ratio (Querkontraktionszahl) v: " << v << "\n";
	ss << "  bulk modulus (Kompressionsmodul): " << kappa << "\n";
	ss << "  Elasticity Tensor = " << elastTensorFunct << "\n";
	m_materialConfiguration = ss.str();
	UG_LOG("\n" << m_materialConfiguration << "\n");
}

}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#endif /* HOOKE_IMPL_H_ */
