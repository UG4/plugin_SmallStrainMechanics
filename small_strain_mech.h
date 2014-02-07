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
* This class assembles the equation of the linear elasticity problem:
*
*  - div(sigma) = f - rho * \frac{\partial^2 u}{\partial^2 t}
*  sigma = C : epsilon
*  epsilon = 1/2 (nabla_u + nabla_u^T)
*
*  As common in the linear theory we identify the deformed configuration with
*  the undeformed configuration. Therefore we do not have to distinguish between
*  different stress measures.
*  Here we use a pure displacement-ansatz to solve the coupled system above.
*
*  \TODO: in case of linear elasticity one can implement the voigt notation (exploiting
*  	  symmetry of sigma and epsilon)
*
*  For plastic material behavior we suppose, that the linearized strain tensor could be
*  decomposed additively:
*
*  eps = eps_e + eps_p.
*
*  The plastic behavior is described by a flow-condition and a flow-rule for the plastic
*  evolution (\frac{\partial eps_p){\partial t} = ...). At present a flow-condition of
*  von-Mises-type and an associative flow-rule are implemented. To treat the plasticity
*  we use the well-established return-mapping-algorithm. This algorithm is valid for the
*  3d-case and the plane strain-case, but not for the plane stress-case!
*  The plastic material behavior should be enabled by the script-call:
*  "elemDisc:use_elastoplast_mat_behavior(true)"
*  before the "elemDisc:set_hooke_elasticity_tensor"-call!
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
		void add_material_law(SmartPtr<IMaterialLaw<TDomain> > spMatLaw)
		{ m_spMatLaw = spMatLaw;}

	///	set an output writer
		void set_output_writer(SmartPtr<MechOutputWriter<TDomain> > spOutWriter){
			m_spOutWriter = spOutWriter;
			m_bOutWriter = true;
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


	///	add mass jacobian
		void add_mass_jacobian(const bool bMassJac){ m_bAddMassJac = bMassJac;}

		/*std::string hardening_config_string() const
		{
			std::stringstream ss;
			switch (m_hardening)
			{
				case 0: return "perfect hardening";
				case 1: ss << "linear hardening, " << matConsts.config_string(); return ss.str();
				case 2:
					ss << "exponential hardening: " << " max hard iter = " << m_MaxHardIter
					   << " hardening accuary = " << m_HardAccuracy << " "
					   << matConsts.config_string();
					return ss.str();
			}
			return "unknown hardening behaviour";
		}*/

	///	sets the quad order
		void set_quad_order(const size_t order) {m_quadOrder = order; m_bQuadOrderUserDef = true;}

	///	initialize state/"inner" variables
		void init_state_variables(const size_t order);

	///	get stress eigenvalues at coord
	/*	void stress_eigenvalues_at(const number CoordX, const number CoordY, const number CoordZ)
		{
			m_stressEV = true;
			m_eigCoordX = CoordX; m_eigCoordY = CoordY; m_eigCoordZ = CoordZ;
		}

	///	get normal stresses at coord
		void normal_stresses_at(const number CoordX, const number CoordY, const number CoordZ)
		{
			m_normalStress = true;
			m_normStressCoordX = CoordX; m_normStressCoordY = CoordY; m_normStressCoordZ = CoordZ;
		}

	/// provides the stress tensor sigma and the linearized strain tensor epsilon
		void normal_stress_strain_loc(LocalVector& devSigma, LocalVector& sigma,
				LocalVector& eps, TBaseElem* elem, const LocalVector& u);

	///	closes the gnuplot-file
		void close_gnuplot_file();
*/

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


	///	returns config information of convergence check and preconditioner
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

			//ss << " Use Plasticity is " << (m_bUsePlasticity ? "ON" : "OFF") << "\n";
			//ss << " Hardening: " << ConfigShift(hardening_config_string()) << "\n";
			//ss << " TangentAccuracy = " << m_tangentAccur << "\n";
			//ss << " Elasticity Configuration: " << ConfigShift(m_materialConfiguration) << "\n";
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

	///	get trace and deviatoric part of a matrix
	//	number MatDeviatorTrace(const MathMatrix<dim, dim>& mat, MathMatrix<dim, dim>& dev);

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

	/// add mass jacobian
		bool m_bAddMassJac;

	///	compute stress-eigenvalues / normalStress at specific corner
		/*bool m_stressEV;
		bool m_normalStress;

	/// coords where to get stress-eigenvales/normalStress
		number m_eigCoordX, m_eigCoordY, m_eigCoordZ;
		number m_normStressCoordX, m_normStressCoordY, m_normStressCoordZ;

	///	checks if ip values are already written, e.g. in order to write stress-eigenvalues
		bool m_bIP_values_written;
*/
	///	maximal value describing the factor used to project back to the elastic domain
	///	(in case of plastic behavior)
		number m_max_gamma;

	///	output-file
	//	FILE *m_outFile;

		std::string m_materialConfiguration;
};

// end group small_strain_mechanics
/// \}

} //end of namespace SmallStrainMechanics
} //end of namespace ug

#endif /* SMALL_STRAIN_MECH_H_ */
