/*
 * small_strain_mechanics_tools_impl.h
 *
 *  Created on: 16.05.2012
 *      Author: raphaelprohl, Andreas Vogel
 */

#include "common/util/string_table_stream.h"

#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

using namespace std;

namespace ug {
namespace SmallStrainMechanics{

/// \addtogroup small_strain_mechanics
/// \{

template<typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
print_mat_constants(const number lambda, const number mu,
		const number E, const number v)
{
	UG_LOG("\n Material parameter: \n");
	UG_LOG( "1. Lame constant lambda: " << lambda << "\n");
	UG_LOG( "2. Lame constant mue (sometimes 'G', shear modulus): " << mu << "\n");
	UG_LOG( "This setting equals: \n");
	UG_LOG( "young modulus E: " << E << "\n");
	UG_LOG( "poisson ratio v: " << v << "\n \n");
}

template<typename TDomain>
number
SmallStrainMechanicsElemDisc<TDomain>::
MatDeviatorTrace(const MathMatrix<dim, dim>& mat, MathMatrix<dim, dim>& dev)
{
	number trace = Trace(mat);

	//	compute the deviatoric part of mat
	for (size_t i = 0; i < (size_t) dim; ++i){
		for (size_t j = 0; j < (size_t) dim; ++j){
			dev[i][j] = mat[i][j];
		}
		dev[i][i] -= 1.0 / dim * trace;
	}

	return trace;
}

template<typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
set_elasticity_tensor_orthotropic(
const number C11, const number C12, const number C13,
		const number C22, const number C23,
					const number C33,
							const number C44,
									const number C55,
											const number C66 )
{
	static const size_t dim = TDomain::dim;

	UG_ASSERT( dim==3, "Orthotrope Tensor only for 3 dimensions" );

	//  setze alle wWerte auf 0
	for (size_t i = 0; i < (size_t) dim; ++i){
		for (size_t j = 0; j < (size_t) dim; ++j){
			for (size_t k = 0; k < (size_t) dim; ++k){
				for (size_t l = 0; l < (size_t) dim; ++l)
				{
					m_ElastTensorFunct[i][j][k][l] = 0.0;
				} //end (l)
			} //end (k)
		} //end (j)
	} //end (i)

	// Tensor mit Werte fuellen
	//                 i  j  k  l
	m_ElastTensorFunct[0][0][0][0] = C11;

	m_ElastTensorFunct[0][0][1][1] = C12;
	m_ElastTensorFunct[1][1][0][0] = C12; // = C21

	m_ElastTensorFunct[0][0][2][2] = C13;
	m_ElastTensorFunct[2][2][0][0] = C13; // = C31

	m_ElastTensorFunct[1][1][1][1] = C22;

	m_ElastTensorFunct[1][1][2][2] = C23;
	m_ElastTensorFunct[2][2][1][1] = C23; // = C32

	m_ElastTensorFunct[2][2][2][2] = C33;

	m_ElastTensorFunct[1][2][1][2] = C44;
	m_ElastTensorFunct[1][2][2][1] = C44;

	m_ElastTensorFunct[2][1][1][2] = C44;
	m_ElastTensorFunct[2][1][2][1] = C44;

	m_ElastTensorFunct[2][0][2][0] = C55;
	m_ElastTensorFunct[0][2][2][0] = C55;
	m_ElastTensorFunct[2][0][0][2] = C55;
	m_ElastTensorFunct[0][2][0][2] = C55;

	m_ElastTensorFunct[0][1][0][1] = C66;
	m_ElastTensorFunct[1][0][0][1] = C66;
	m_ElastTensorFunct[0][1][1][0] = C66;
	m_ElastTensorFunct[1][0][1][0] = C66;

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


template<typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
set_hooke_elasticity_tensor_E_nu(const number E, const number nu)
{
	number lambda = (E*nu) / ((1+nu)*(1-2*nu));
	number mu = E/(2*(1+nu));
	set_hooke_elasticity_tensor(lambda, mu);

}

template<typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
set_hooke_elasticity_tensor(const number lambda, const number mu)
{
	//	sets the 'Hooke'-elasticity tensor for isotropic materials
	//	in 2D this tensor formulation corresponds to the
	//	plane strain assumption for Hookes`s law
	static const int dim = TDomain::dim;

	matConsts.mu = mu;
	matConsts.kappa = lambda + 2.0/3.0 * mu;

	//  filling the constant elasticity tensor
	for (size_t i = 0; i < (size_t) dim; ++i){
		for (size_t j = 0; j < (size_t) dim; ++j){
			for (size_t k = 0; k < (size_t) dim; ++k){
				for (size_t l = 0; l < (size_t) dim; ++l)
				{
					m_ElastTensorFunct[i][j][k][l] = 0.0;

					if ((i == j) && (k == l)) {
						m_ElastTensorFunct[i][j][k][l] += lambda;
					}

					if ((i == k) && (j == l)) {
						m_ElastTensorFunct[i][j][k][l] +=  mu;
					}

					if ((i == l) && (j == k)) {
						m_ElastTensorFunct[i][j][k][l] +=  mu;
					}
				} //end (l)
			} //end (k)
		} //end (j)
	} //end (i)

	number E = mu * (3.0 * lambda + 2.0 * mu)/(lambda + mu);
	number v = 0.5 * lambda/(lambda + mu);
	number K = (3*lambda + 2*mu)/2;

	std::stringstream ss;
	ss << "Hooke Elasticity Tensor: \n";
	ss << "  Lame constant lambda: " << lambda << "\n";
	ss << "  Lame constant mue (sometimes 'G', shear modulus): " << mu << "\n";
	ss << " This setting equals: \n";
	ss << "  young modulus E: " << E << "\n";
	ss << "  poisson ratio v: " << v << "\n";
	ss << "  bulk modulus (Kompressionsmodul) " << K << "\n";
	ss << "  Elasticity Tensor = " << m_ElastTensorFunct << "\n";
	m_materialConfiguration = ss.str();
	UG_LOG("\n" << m_materialConfiguration << "\n");
}

//	"DisplacementGradient": computes the displacement-gradient at a ip
template<typename TDomain>
template<typename TFEGeom>
void
SmallStrainMechanicsElemDisc<TDomain>::
DisplacementGradient(MathMatrix<dim, dim>& GradU, const TFEGeom& geo,
		const LocalVector& u, const size_t ip)
{
	//	loop shape-functions at one integration point ip in order
	//	to compute local_grad(ip,i): \frac{\partial N_i}{\eps_ip}
	//	and global_grad(ip,i): \frac{\partial N_i}{\X_ip}

	for (size_t i = 0; i < (size_t) dim; ++i) {
		for (size_t J = 0; J < (size_t) dim; ++J)
		{
			GradU[i][J] = 0.0;
			for (size_t a = 0; a < geo.num_sh(); ++a) // loop shape-fcts
			{
				//	compute GradU: displacementGradient
				GradU[i][J] += geo.global_grad(ip, a)[J] * u(i, a);
			}
		}
	}
}

/*
 *  for more info about the hardening-behavior see Simo/Hughes 98 p.121,122 for example!
 */
//	"Hardening": nonlinear hardening law
template<typename TDomain>
inline number
SmallStrainMechanicsElemDisc<TDomain>::
Hardening(const number alpha)
{
	return (matConsts.K_0 + matConsts.Hard * alpha +
			(matConsts.K_inf - matConsts.K_0) * (1.0 - exp(-matConsts.omega * alpha)));
}

//	"Hardening_d": derivative of nonlinear hardening law "Hardening"
template<typename TDomain>
inline number
SmallStrainMechanicsElemDisc<TDomain>::
Hardening_d(const number alpha)
{
	return (matConsts.Hard + (matConsts.K_inf - matConsts.K_0) * matConsts.omega
			* exp(-matConsts.omega * alpha));
}
template<typename TDomain>
inline number
SmallStrainMechanicsElemDisc<TDomain>::
PerfectPlasticity(const number flowcondtrial)
{
	number gamma = flowcondtrial / (2.0 * matConsts.mu);
	return gamma;
}

template<typename TDomain>
inline number
SmallStrainMechanicsElemDisc<TDomain>::
LinearHardening(const number flowcondtrial)
{
	number gamma = flowcondtrial / (2.0 * (matConsts.mu + matConsts.Hard/3.0));
	return gamma;
}

template<typename TDomain>
number
SmallStrainMechanicsElemDisc<TDomain>::
ExponentialHardening(const number strialnorm,
		const number alpha, const number mu)
{
	number gamma = 0.0;

	for (size_t i = 0; i < m_MaxHardIter; ++i) {
		//  f_cap has to be equal to 0, due to the consistency condition
		number f_cap = strialnorm - 2.0 * gamma * mu - sqrt(2.0 / 3.0) * Hardening(
				alpha + sqrt(2.0 / 3.0) * gamma);

		//abs(f_cap) < HardAccuracy * strialnorm vs. - f_cap > - HardAccuracy * strialnorm
		if (f_cap < m_HardAccuracy * strialnorm){ break;}

		gamma += f_cap / (2.0 * mu * (1.0 + Hardening_d(
				alpha + sqrt(2.0 / 3.0) * gamma) / (3.0 * mu)));
	}

	return gamma;
}

template<typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
Update(const MathMatrix<dim, dim>& GradU,
		const MathMatrix<dim, dim>& eps_p_old_t,
		MathMatrix<dim, dim>& eps_p_new,
		number& alpha)
{
	//	compute linearized strain tensor (eps)
	MathMatrix<dim, dim> eps;
	for(size_t i = 0; i < (size_t) dim; ++i){
		for(size_t j = 0; j < (size_t) dim; ++j){
			eps[i][j] = 0.5 * (GradU[i][j] + GradU[j][i]);
		}
	}
	MathMatrix<dim, dim> dev_eps;
	MatDeviatorTrace(eps, dev_eps);

	MathMatrix<dim, dim> strial;
	for(size_t i = 0; i < (size_t) dim; ++i){
		for(size_t j = 0; j < (size_t) dim; ++j){
			strial[i][j] = 2.0 * matConsts.mu * (dev_eps[i][j] - eps_p_old_t[i][j]);
		}
	}

	number strialnorm = MatFrobeniusNorm(strial);
	number flowcondtrial = strialnorm - sqrt(2.0 / 3.0) * Hardening(alpha);

	if (flowcondtrial <= 0){
		eps_p_new = eps_p_old_t; return;
	}

	number gamma = 0.0;

	UG_ASSERT(strialnorm > 0.0, "norm of strial needs to be > 0.0");

	/*	computation of gamma (plastic corrector/multiplicator) */
	switch (m_hardening)
	{
		case 0: gamma = PerfectPlasticity(flowcondtrial); break;
		case 1: gamma = LinearHardening(flowcondtrial); break;
		case 2: gamma = ExponentialHardening(strialnorm, alpha, matConsts.mu); break;
		default:
			UG_THROW(m_hardening << " in 'Update' is not a valid hardening behavior! \n");
	}

	MathMatrix<dim, dim> normaltrial;
	MatScale(normaltrial, 1.0 / strialnorm, strial);

	UG_ASSERT(gamma > 0.0, "gamma needs to be > 0.0");

	MathMatrix<dim, dim> d_eps_p;
	MatScale(d_eps_p, gamma, normaltrial);
	MatAdd(eps_p_new, eps_p_old_t, d_eps_p);

	alpha += sqrt(2.0 / 3.0) * gamma;

	if (abs(gamma) > abs(m_max_gamma))
		m_max_gamma = gamma;
}

template<typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
Flowrule(const MathMatrix<dim, dim>& GradU,
		const MathMatrix<dim, dim>& eps_p_old_t,
		const number alpha,
		MathMatrix<dim, dim>& strial,
		number& gamma,
		MathMatrix<dim, dim>& normal,
		MathMatrix<dim, dim>& eps,
		MathMatrix<dim, dim>& eps_p_new)
{
	//////////////////////////
	//  TRIAL ELASTIC STEP
	//////////////////////////

	MathMatrix<dim, dim> eps_trial;

	//	compute linearized strain tensor (eps)
	for(size_t i = 0; i < (size_t) dim; ++i){
		for(size_t j = 0; j < (size_t) dim; ++j)
		{
			eps[i][j] = 0.5 * (GradU[i][j] + GradU[j][i]);
			//	get eps_trial := eps_{n+1} - eps^p_{n}
			eps_trial[i][j] = eps[i][j] - eps_p_old_t[i][j];
		}
	}
	MathMatrix<dim, dim> dev_eps_trial;
	MatDeviatorTrace(eps_trial, dev_eps_trial);

	//	compute trial strain deviator
	MatScale(strial, 2.0 * matConsts.mu, dev_eps_trial);

	number strialnorm = MatFrobeniusNorm(strial);
	number flowcondtrial = strialnorm - sqrt(2.0 / 3.0) * Hardening(alpha);

	//////////////////////////////
	//  CHECKING YIELD-CONDITION:
	//////////////////////////////

	if (flowcondtrial <= 0){
		eps_p_new = eps_p_old_t; gamma = 0.0; normal = 0.0;
		return;
	}

	////////////////////////////////////
	//  RETURN-MAPPING (corrector-step)
	////////////////////////////////////

	UG_ASSERT(strialnorm > 0.0, "norm of strial needs to be > 0.0");

	/*	computation of gamma (plastic corrector/multiplicator)
	 * 	accordingly to Simo/Hughes 98 p.121/122 */
	switch (m_hardening)
	{
		case 0: gamma = PerfectPlasticity(flowcondtrial); break;
		case 1: gamma = LinearHardening(flowcondtrial); break;
		case 2: gamma = ExponentialHardening(strialnorm, alpha, matConsts.mu); break;
		default:
			UG_THROW(m_hardening << " in 'Flowrule' is not a valid hardening behavior! \n");
	}

	UG_ASSERT(gamma > 0.0, "gamma needs to be > 0.0");

	MatScale(normal, 1.0 / strialnorm, strial);

	MathMatrix<dim, dim> d_eps_p;
	MatScale(d_eps_p, gamma, normal);
	MatAdd(eps_p_new, eps_p_old_t, d_eps_p);
}

template<typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
ConstLaw(const MathMatrix<dim, dim>& eps, const MathMatrix<dim, dim>& strial,
		const number& gamma, const MathMatrix<dim, dim>& normal,
		MathMatrix<dim, dim>& T)
{
	//	constitutive law taken from Simo/Hughes 98 'Computational Inelasticity' p. 124
	number trEps = Trace(eps);

	//	compute sigma = kappa * tr[eps] * id + strial - 2 * mu * gamma * normal
	for(size_t i = 0; i < (size_t) dim; ++i){
		for(size_t j = 0; j < (size_t) dim; ++j){
			T[i][j] = strial[i][j] - 2.0 * matConsts.mu * gamma * normal[i][j];
		}
		T[i][i] += matConsts.kappa * trEps;
	}

}


template<typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
StressTensor(MathMatrix<dim, dim>& sigma, const MathMatrix<dim, dim>& GradU,
		const MathMatrix<dim, dim>& eps_p_old_t, const number alpha)
{
	MathMatrix<dim, dim> strial, normal, eps, eps_p_new;
	number gamma;

	Flowrule(GradU, eps_p_old_t, alpha, strial, gamma, normal, eps, eps_p_new);
	ConstLaw(eps, strial, gamma, normal, sigma);
}


//	elasticity tensor by means of numerical differentiation
template<typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
TangentNumApprox(MathTensor4<dim, dim, dim, dim>& C, MathMatrix<dim, dim> GradU,
		const MathMatrix<dim, dim>& eps_p_old_t, const number alpha)
{
	MathMatrix<dim, dim> T, TT;

	//	for formulation in reference config
	StressTensor(T, GradU, eps_p_old_t, alpha);

	for (size_t i = 0; i < (size_t) dim; ++i){
		for (size_t j = 0; j < (size_t) dim; ++j)
		{
			GradU[i][j] += m_tangentAccur;

			StressTensor(TT, GradU, eps_p_old_t, alpha);

			for (size_t k = 0; k < (size_t) dim; ++k){
				for (size_t l = 0; l < (size_t) dim; ++l){
					C[i][j][k][l] = 1.0/ m_tangentAccur * (TT[k][l] - T[k][l]);
				}
			}

			GradU[i][j] -= m_tangentAccur;
		}
	}

}

/* this function computes the contraction of a 4th order tensor by a second order tensor */
template<typename TDomain>
void SmallStrainMechanicsElemDisc<TDomain>::
TensContract4(MathMatrix<dim, dim>& m_out, const MathTensor4<dim, dim, dim, dim>& tens4,
		const MathMatrix<dim, dim>& tens2)
{
	for(size_t i = 0; i < (size_t) dim; ++i){
		for(size_t j = 0; j < (size_t) dim; ++j)
		{
			m_out[i][j] = 0.0;

			for(size_t k = 0; k < (size_t) dim; ++k){
				for(size_t l = 0; l < (size_t) dim; ++l)
				{
					m_out[i][j] += tens4[i][j][k][l] * tens2[k][l];
				}
			}
		}
	}
}

// end group small_strain_mechanics
/// \}

} //end of namespace SmallStrainMechanics
} //end of namespace ug
