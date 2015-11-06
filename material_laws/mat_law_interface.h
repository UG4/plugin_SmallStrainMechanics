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


		virtual SmartPtr<MathMatrix<dim, dim> > inelastic_strain_tensor(const size_t ip)
		{
			MathMatrix<dim, dim> inelastStrain; inelastStrain = 0;
			SmartPtr<MathMatrix<dim, dim> > spInelasticStrain (new MathMatrix<dim, dim>(inelastStrain));
			return spInelasticStrain;
		}
		virtual number hardening_parameter(const size_t ip){return 0.0;}
		virtual number plastic_multiplier(const size_t ip, const MathMatrix<dim, dim>& GradU)
			{return 0.0;}

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
