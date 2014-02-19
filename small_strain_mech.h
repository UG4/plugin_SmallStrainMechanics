/*
 * small_strain_mech.h
 *
 *  Created on: 16.05.2012
 *      Author: raphaelprohl, Andreas Vogel
 */

#ifndef SMALL_STRAIN_MECH_H_
#define SMALL_STRAIN_MECH_H_

#include <stdio.h>
#include <string>

#include "bridge/bridge.h"

#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/local_finite_element/lagrange/lagrange.h"
#include "lib_disc/quadrature/gauss/gauss_quad.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "material_laws/mat_law_interface.h"
#include "output_writer/mech_output_writer.h"

namespace ug{
namespace SmallStrainMechanics{

/// \addtogroup small_strain_mechanics
/// \{
/// Element Discretization for Linear Elasticity and Elasto-Plasticity problems (restricted to small deformations)
/**
* This class assembles the equation of the static linear elasticity problem:
*
*  (1)		- div(sigma) = f
*  (2)		sigma = C : epsilon
*  (3)		epsilon = 1/2 (nabla_u + nabla_u^T)
*
*		+ boundary conditions.
*
*  As common in the linear theory we identify the deformed configuration with
*  the undeformed configuration. Therefore we do not have to distinguish between
*  different stress measures.
*  Here we use a pure displacement-ansatz to solve the coupled system above,
*  i.e. the kinematic equation for the linearized strain tensor epsilon (3) is inserted in
*  the material law (2) and the resulting stress tensor sigma is introduced in the momentum balance
*  equation (1). Therein, the computation of the strain and the stress tensor is performed at the
*  integration points! Following this approach, the only remaining unknown function
*  is the displacement field u!
*
*  \TODO: in case of linear elasticity one can implement the voigt notation (exploiting
*  	  symmetry of sigma and epsilon)
*
*  In order to compute the eigenfrequencies of a system by means of the eigenvalue problem
*  the second time derivatives of u
*
*  		- rho * \frac{\partial^2 u}{\partial^2 t}
*
*  are considered in the momentum equation (1) and deliver a contribution to the mass matrix.
*
* References:
* <ul>
* <li> M. Rupp. Berechnung der Resonanzschwingungen einer Gitarrendecke.
* <li> (Diplomarbeit, 2009, Universität Heidelberg)
* <ul>
*/

template <typename TDomain>
class SmallStrainMechanicsElemDisc
	: public IElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

	///	own type
		typedef SmallStrainMechanicsElemDisc<TDomain> this_type;

	///	base element type of associated domain
		typedef typename domain_traits<TDomain::dim>::geometric_base_object TBaseElem;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	public:
	/// constructor
		SmallStrainMechanicsElemDisc(const char* functions, const char* subsets);

	///	destructor
		virtual	~SmallStrainMechanicsElemDisc();

	///	adds a material law
		void set_material_law(SmartPtr<IMaterialLaw<TDomain> > spMatLaw)
		{ m_spMatLaw = spMatLaw;}

	///	set an output writer
		void set_output_writer(SmartPtr<MechOutputWriter<TDomain> > spOutWriter){
			m_spOutWriter = spOutWriter; m_bOutWriter = true;
			m_spOutWriter->preprocess();
		}

	///	set volume forces
		void set_volume_forces(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
		void set_volume_forces(number vel);
		void set_volume_forces(number vel_x, number vel_y);
		void set_volume_forces(number vel_x, number vel_y, number vel_z);
	#ifdef UG_FOR_LUA
			void set_volume_forces(const char* fctName);
	#endif

		/**
		 * This method sets the pressure term. A zero value is assumed as default.
		 */
		///	\{
		void set_pressure(SmartPtr<CplUserData<number, dim> > user);
		void set_pressure(number val);
#ifdef UG_FOR_LUA
		void set_pressure(const char* fctName);
#endif
		///	\}

	///	sets the quad order
		void set_quad_order(const size_t order) {m_quadOrder = order; m_bQuadOrderUserDef = true;}

	///	initialize state/"inner" variables
		void init_state_variables(const size_t order);

	/// computing contact forces elementwise by averaging over all integration points
		template<typename TElem>
		void contact_forces_elem_ips_avg(LocalVector& locForce, GeometricObject* side, TElem* elem,
				const MathVector<dim> sideCoPos[], int numElemCorners, const LocalVector& locU,
				std::vector<DoFIndex> vActiveSetLoc);

	/// computing contact forces elementwise at every element midpoint
		template <typename TSide, typename TElem>
		void contact_forces_elem_midpoint(LocalVector& locForce, TSide* side, TElem* elem,
				const MathVector<dim> sideCoPos[], const LocalVector& locU,
				std::vector<DoFIndex> vActiveSetLoc);


	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
		//	check number
			if(vLfeID.size() != (size_t)dim)
				UG_THROW("SmallStrainMechanics: Needs exactly "<<dim<<" functions.");

			//	check & set order
			for(size_t i = 0; i < vLfeID.size(); ++i)
				if(vLfeID[i].order() < 1)
					UG_THROW("SmallStrainMechanics: Adaptive order not implemented.");

		//	remember lfeID;
			m_lfeID = vLfeID[0];

		//	set order
			m_order = vLfeID[0].order();

		//	update assemble functions
			set_assemble_funcs();
		}


	///	returns config information of small strain mechanics ElemDisc and material law
		std::string config_string() const
		{
			std::stringstream ss;
			ss << "SmallStrainMechanics " << dim << "d ";
			if(dim == 2)
				ss << " [Plain Strain / Ebener Verzerrungszustand]";
			ss << "\n";
			ss << " order of disc scheme = " << m_order << "\n";
			ss << " shape function set " << m_lfeID << "\n";
			if(m_bQuadOrderUserDef)
				ss << " User Defined Quad Order = " << m_quadOrder << "\n";

			ss << " Material Configuration: " << ConfigShift(m_spMatLaw->m_materialConfiguration) << "\n";
			return ss.str();
		}

	private:
		size_t num_fct() const {return dim;}

		///	assemble methods
		template<typename TElem, typename TFEGeom>
		void prep_timestep_elem(const number time, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		template<typename TElem, typename TFEGeom>
		void prep_elem(const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void fsh_elem_loop();

		template<typename TElem, typename TFEGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void add_rhs_elem(LocalVector& d, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void fsh_timestep_elem(const number time, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

	protected:
		template <typename TElem, typename TFEGeom>
		void lin_def_pressure(const LocalVector& u,
				              std::vector<std::vector<number> > vvvLinDef[],
				              const size_t nip);

		template <typename TElem, typename TFEGeom>
		void lin_def_volume_forces(const LocalVector& u,
				                  std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
				                  const size_t nip);

		///	computes the displacements (and derivatives)
		template <typename TElem, typename TFEGeom>
		void ex_displacement(const LocalVector& u,
				          	 const MathVector<dim> vGlobIP[],
				             const MathVector<TFEGeom::Type::dim> vLocIP[],
				             const size_t nip,
				             MathVector<dim> vValue[],
				             bool bDeriv,
				             std::vector<std::vector<MathVector<dim> > > vvvDeriv[]);

	private:
	///	updates the geometry for a given element
		void update_geo_elem(TBaseElem* elem, DimFEGeometry<dim>& geo);

	///	prints material constants
		void print_mat_constants(const number lambda,
				const number mu, const number E, const number v);

	private:
		///	sets the requested assembling routines
		void set_assemble_funcs();

		void register_all_fe_funcs(int order, int quadOrder);
		template <typename TElem, typename TFEGeom>
		void register_fe_func();

	private:
	///	material law
		SmartPtr<IMaterialLaw<TDomain> > m_spMatLaw;

	///	elasticity tensor
		SmartPtr<MathTensor4<dim,dim,dim,dim> > m_spElastTensor;

	///	output writer
		SmartPtr<MechOutputWriter<TDomain> > m_spOutWriter;
		bool m_bOutWriter;
		bool m_bMatLawPassedToOutWriter;

	///	current order of disc scheme
		int m_order;

	///	current shape function set
		LFEID m_lfeID;

	///	current integration order
		bool m_bQuadOrderUserDef;
		int m_quadOrder;

	///	data import for volume forces
		DataImport<MathVector<dim>, dim > m_imVolForce;

	///	Data import for the reaction term
		DataImport<number, dim> m_imPressure;
};

// end group small_strain_mechanics
/// \}

} //end of namespace SmallStrainMechanics
} //end of namespace ug

#endif /* SMALL_STRAIN_MECH_H_ */
