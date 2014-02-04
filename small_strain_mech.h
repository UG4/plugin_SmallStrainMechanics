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
		void set_hooke_elasticity_tensor_E_nu(const number E,const number nu);

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


	///	use elastoplastic material behavior
		void use_elastoplast_mat_behavior(const bool elastplast){ m_bUsePlasticity = elastplast;}

		std::string hardening_config_string() const
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
		}
	///	set hardening behavior
		void set_hardening_behavior(int hard)
		{
			switch (hard)
			{
				//	perfect hardening
				case 0: m_hardening = 0; break;

				//	linear hardening
				case 1: m_hardening = 1; matConsts.Hard = 129.24; break;

				//	exponential hardening
				case 2: m_hardening = 2; matConsts.Hard = 129.24;
						matConsts.omega = 16.93; m_MaxHardIter = 100;
						m_HardAccuracy = 1e-10; break;

				default: UG_THROW(hard << " is not a valid hardening behavior! "
						"Choose 0 (perfect), 1 (linear) or 2 (exponential) ! \n");
			}

			matConsts.K_0 = 450.00; //= 450e6 Pa
			matConsts.K_inf = 715.00;
		}

	///	use numerical approximation the tangent; tanAccur: tangent-accuracy
		void use_approx_tangent(const number tanAccur){m_tangentAccur = tanAccur;}

	///	sets the quad order
		void set_quad_order(const size_t order) {m_quadOrder = order; m_bQuadOrderUserDef = true;}

	///	initialize state/"inner" variables
		void init_state_variables(const size_t order);

	///	get stress eigenvalues at coord
		void stress_eigenvalues_at(const number CoordX, const number CoordY, const number CoordZ)
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

			ss << " Use Plasticity is " << (m_bUsePlasticity ? "ON" : "OFF") << "\n";
			ss << " Hardening: " << ConfigShift(hardening_config_string()) << "\n";
			ss << " TangentAccuracy = " << m_tangentAccur << "\n";
			ss << " Elasticity Configuration: " << ConfigShift(m_materialConfiguration) << "\n";
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


	protected:
		virtual void approximation_space_changed()
		{
			clear_attachments();
			attach_attachments();
		}


	private:
	///	updates the geometry for a given element
		void update_geo_elem(TBaseElem* elem, DimFEGeometry<dim>& geo);

	///	prints material constants
		void print_mat_constants(const number lambda,
				const number mu, const number E, const number v);

	///	get trace and deviatoric part of a matrix
		number MatDeviatorTrace(const MathMatrix<dim, dim>& mat, MathMatrix<dim, dim>& dev);

		template<typename TFEGeom>
		void DisplacementGradient(MathMatrix<dim, dim>& GradU, const TFEGeom& geo,
				const LocalVector& u, const size_t ip);

		number Hardening(const number alpha);

		number Hardening_d(const number alpha);

		number PerfectPlasticity(const number flowcondtrial);

		number LinearHardening(const number flowcondtrial);

		number ExponentialHardening(const number strialnorm,
				const number alpha, const number mu);

		void Update(const MathMatrix<dim, dim>& GradU, const MathMatrix<dim, dim>& eps_p_old_t,
				MathMatrix<dim, dim>& eps_p_new, number& alpha);

		void Flowrule(const MathMatrix<dim, dim>& GradU, const MathMatrix<dim, dim>& eps_p_old_t,
				const number alpha, MathMatrix<dim, dim>& strial, number& gamma,
				MathMatrix<dim, dim>& normal, MathMatrix<dim, dim>& eps,
				MathMatrix<dim, dim>& eps_p_new);

		void ConstLaw(const MathMatrix<dim, dim>& eps, const MathMatrix<dim, dim>& strial,
				const number& gamma, const MathMatrix<dim, dim>& normal, MathMatrix<dim, dim>& T);

		void StressTensor(MathMatrix<dim, dim>& sigma, const MathMatrix<dim, dim>& GradU,
				const MathMatrix<dim, dim>& eps_p_old_t, const number alpha);

		void TangentNumApprox(MathTensor4<dim, dim, dim, dim>& C, MathMatrix<dim, dim> GradU,
				const MathMatrix<dim, dim>& eps_p_old_t, const number alpha);

		void TensContract4(MathMatrix<dim, dim>& m_out,
				const MathTensor4<dim, dim, dim, dim>& tens4, const MathMatrix<dim, dim>& tens2);


	private:
	///	use this method to make sure that all required attachments are attached
	/**	This method won't be necessary if we attach m_aElemData during initialization.*/
		void attach_attachments()
		{
			typename TDomain::grid_type& grid = *this->domain()->grid();
	//todo:	Move this to an initialization routine. Best would probably to make
	//		set_domain(...) of the base class virtual and to attach m_aElemData there.
			grid.template attach_to<TBaseElem>(m_aElemData);
			m_aaElemData.access(grid, m_aElemData);
		}

		void clear_attachments()
		{
			typename TDomain::grid_type& grid = *this->domain()->grid();
			if(grid.template has_attachment<TBaseElem>(m_aElemData)){
				grid.template detach_from<TBaseElem>(m_aElemData);
				m_aaElemData.invalidate();
			}
		}

		///	sets the requested assembling routines
		void set_assemble_funcs();

		void register_all_fe_funcs(int order, int quadOrder);
		template <typename TElem, typename TFEGeom>
		void register_fe_func();

	private:
	///	current order of disc scheme
		int m_order;

	///	current shape function set
		LFEID m_lfeID;

	///	current integration order
		bool m_bQuadOrderUserDef;
		int m_quadOrder;

	/// elasticity tensor
		MathTensor4<dim, dim, dim, dim> m_ElastTensorFunct;

	///	data import for volume forces
		DataImport<MathVector<dim>, dim > m_imVolForce;

	///	Data import for the reaction term
		DataImport<number, dim> m_imPressure;

	/// elastoplastic material behavior
		bool m_bUsePlasticity;

	/// hardening behavior
		int m_hardening;
		size_t m_MaxHardIter;
		number m_HardAccuracy;

		struct MaterialConstants{
			number mu;		//	shear modulus
			number kappa;	//	bulk modulus
			number K_0; 	//	initial flow-stress:
			number K_inf; 	//	residual flow-stress (saturation stress)
			number Hard;	//	linear hardening modulus
			number omega;	//	hardening exponent

			std::string config_string() const
			{
				std::stringstream ss;
				ss << "MaterialConstant:\n";
				if(mu != 0.0) ss << " shear modulus mu = " << mu << "\n";
				if(kappa != 0.0) ss << " bulk modulus kappa = " << kappa << "\n";
				if(K_inf != 0.0) ss << " residual flow-stress (saturation stress) K_inf = " << K_inf << "\n";
				if(Hard != 0.0)  ss << " linear hardening modulus = " << Hard << "\n";
				if(omega != 0.0) ss << " hardening exponent omega = " << omega;
				return ss.str();
			}
		}matConsts;

	/// tangent accuracy
		number m_tangentAccur;

	///	compute stress-eigenvalues / normalStress at specific corner
		bool m_stressEV;
		bool m_normalStress;

	/// coords where to get stress-eigenvales/normalStress
		number m_eigCoordX, m_eigCoordY, m_eigCoordZ;
		number m_normStressCoordX, m_normStressCoordY, m_normStressCoordZ;

	///	checks if ip values are already written, e.g. in order to write stress-eigenvalues
		bool m_bIP_values_written;

	///	maximal value describing the factor used to project back to the elastic domain
	///	(in case of plastic behavior)
		number m_max_gamma;

	///	output-file
		FILE *m_outFile;

	/// attached ElemData
		struct PlasticUserData {
			MathMatrix<dim, dim> eps_p_old_t; 	//plastic part of strain tensor wrt to last time step
			number alpha;						//hardening variable
		};

		struct ElemData{
			std::vector<PlasticUserData> data; 	//std-vector of PlasticUserData
		};

		ElemData* m_pElemData;

		typedef Attachment<ElemData> AElemData; 	//attachment type: attachment of ElemDatas
		AElemData m_aElemData;						//the instance of the attachment type
		typedef Grid::AttachmentAccessor<TBaseElem, AElemData>	ElemDataAccessor;
		ElemDataAccessor m_aaElemData;

		std::string m_materialConfiguration;
};

// end group small_strain_mechanics
/// \}

} //end of namespace SmallStrainMechanics
} //end of namespace ug

#endif /* SMALL_STRAIN_MECH_H_ */
