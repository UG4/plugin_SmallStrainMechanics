/*
 * hooke_law.h
 *
 *  Created on: 03.02.2014
 *      Author: raphaelprohl
 */

#ifndef HOOKE_LAW_H_
#define HOOKE_LAW_H_

namespace ug{
namespace SmallStrainMechanics{

template <int dim>
class HookeLaw
	: public IMaterialLaw<dim>
{
	public:
	///	constructor
		HookeLaw(){};

	///	Destructor
		~HookeLaw(){};

	protected:
	////////////////////////////
	// INTERFACE-METHODS
	////////////////////////////
		void init();

	///	computes the cauchy stress tensor sigma at an integration point 'ip'
	///	within an element 'elem'
		void stressTensor(MathMatrix<dim,dim>& stressTens, const LocalVector& u, const MathMatrix<dim, dim>& GradU,
				const size_t ip, GeometricObject* elem);

	///	computes the elasticity tensor; commonly denoted by C
		void elasticityTensor(MathTensor4<dim,dim,dim,dim>& elastTens,
				const LocalVector& u, const size_t ip, GeometricObject* elem);

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

	///	returns config information
		/*std::string config_string() const
		{
			std::stringstream ss;
			ss << "HookeLaw " << dim << "d ";
			if(dim == 2)
				ss << " [Plain Strain / Ebener Verzerrungszustand]";
			ss << " Elasticity Configuration: " << ConfigShift(m_materialConfiguration) << "\n";
			return ss.str();
		}*/

	private:
	/// elasticity tensor
		MathTensor4<dim, dim, dim, dim> m_ElastTensorFunct;

		struct MaterialConstants{
			number mu;		//	shear modulus
			number kappa;	//	bulk modulus

			/*std::string config_string() const
			{
				std::stringstream ss;
				ss << "MaterialConstant:\n";
				if(mu != 0.0) ss << " shear modulus mu = " << mu << "\n";
				if(kappa != 0.0) ss << " bulk modulus kappa = " << kappa << "\n";
				return ss.str();
			}*/
		}matConsts;

		//std::string m_materialConfiguration;
};

}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#include "hooke_law_impl.h"

#endif /* HOOKE_LAW_H_ */
