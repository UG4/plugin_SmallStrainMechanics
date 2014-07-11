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

//	NOTE: TDomain is necessary here, due to the attached ElemData, only!
//	-> if an additional interface for materialLaws with internal Vars will be implemented
//	remove TDomain!
template <typename TDomain>
class IMaterialLaw
{
	public:
	///	World dimension
		static const int dim = TDomain::dim;

	///	base element type of associated domain
		typedef typename domain_traits<TDomain::dim>::grid_base_object TBaseElem;

	public:
	///	constructor
		IMaterialLaw(): m_bInit(false), m_bConstElastTens(false){}

	///	destructor
		virtual ~IMaterialLaw(){};

	public:
		///////////////////////////////////////
		//	methods for common material laws
		///////////////////////////////////////
		virtual void init() = 0;

		//	computes a stress tensor at an integration point ip
		virtual void stressTensor(MathMatrix<dim,dim>& stressTens, const size_t ip,
				const MathMatrix<dim,dim>& GradU) = 0;

		//	computes the elasticity tensor at an integration point ip
		virtual SmartPtr<MathTensor4<dim,dim,dim,dim> >
			elasticityTensor(const size_t ip, const MathMatrix<dim, dim>& GradU) = 0;

		//	computes the constant elasticity tensor
		virtual SmartPtr<MathTensor4<dim,dim,dim,dim> > elasticityTensor(){
			SmartPtr<MathTensor4<dim,dim,dim,dim> > spElastTens(new MathTensor4<dim,dim,dim,dim>());
			return spElastTens;
		};

		virtual bool needs_to_add_jac_m(){return true;}

		///////////////////////////////////////
		//	methods for material laws
		//	with internal variables
		///////////////////////////////////////
		virtual void attach_internal_vars(typename TDomain::grid_type& grid){};
		virtual void clear_attachments(typename TDomain::grid_type& grid){};

		virtual void init_internal_vars(TBaseElem* elem, const size_t numIP){};
		virtual void internal_vars(TBaseElem* elem){};
		virtual void update_internal_vars(const size_t ip, const MathMatrix<dim, dim>& GradU){};

		virtual void write_data_to_console(const number t){};

	public:
		inline bool is_initialized(){return m_bInit;}

		inline bool elastTensIsConstant(){return m_bConstElastTens;}

		template <typename TFEGeom>
		void DisplacementGradient(MathMatrix<dim, dim>& GradU, const size_t ip,
				const TFEGeom& geo, const LocalVector& u);

	public:
		std::string m_materialConfiguration;

	protected:
	///	flag indicating, if material law has been initialized
		bool m_bInit;

	/// flag indicating, if elasticity tensor is constant
		bool m_bConstElastTens;

};

}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#include "mat_law_interface_impl.h"

#endif /* MAT_LAW_INTERFACE_H_ */
