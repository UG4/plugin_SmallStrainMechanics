/*
 * mat_law_interface.h
 *
 *  Created on: 03.02.2014
 *      Author: raphaelprohl
 */

#ifndef MAT_LAW_INTERFACE_H_
#define MAT_LAW_INTERFACE_H_

namespace ug{
namespace SmallStrainMechanics{

template <int dim>
class IMaterialLaw
{
	public:
	///	constructor
		IMaterialLaw(){}

	///	destructor
		virtual ~IMaterialLaw(){};

	protected:
		///////////////////////////////////////
		//	methods for common material laws
		///////////////////////////////////////
		virtual void init() = 0;

		//	computing a stress tensor at an integration point ip
		//	of element 'elem'
		virtual void stressTensor(MathMatrix<dim,dim>& stressTens, const LocalVector& u,
				const MathMatrix<dim,dim>& GradU, const size_t ip, GeometricObject* elem) = 0;

		//	computing the elasticity tensor at an integration point ip
		//	of element 'elem'
		virtual void elasticityTensor(MathTensor4<dim,dim,dim,dim>& elastTens,
				const LocalVector& u, const size_t ip, GeometricObject* elem) = 0;

		virtual void print(){};

		///////////////////////////////////////
		//	methods for material laws
		//	with internal variables
		///////////////////////////////////////
		virtual void elem_data(GeometricObject* elem){};
		virtual void update_elem_data(GeometricObject* elem){};

	protected:
		template <typename TFEGeom>
		void DisplacementGradient(MathMatrix<dim, dim>& GradU, const TFEGeom& geo,
				const LocalVector& u, const size_t ip);
};

}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#include "mat_law_interface_impl.h"

#endif /* MAT_LAW_INTERFACE_H_ */
