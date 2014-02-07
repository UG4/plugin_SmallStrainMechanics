/*
 * small_strain_mech.cpp
 *
 *  Created on: 16.05.2012
 *      Author: raphaelprohl, Andreas Vogel
 */

#include "small_strain_mech.h"
#include "small_strain_mech_tools_impl.h"

// for various user data
#include "bindings/lua/lua_user_data.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"

#define PROFILE_SMALL_STRAIN_MECH
#ifdef PROFILE_SMALL_STRAIN_MECH
	#define SMALL_STRAIN_MECH_PROFILE_FUNC()		PROFILE_FUNC_GROUP("Small Strain Mech")
	#define SMALL_STRAIN_MECH_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "Small Strain Mech")
	#define SMALL_STRAIN_MECH_PROFILE_END()		PROFILE_END()
#else
	#define SMALL_STRAIN_MECH_PROFILE_FUNC()
	#define SMALL_STRAIN_MECH_PROFILE_BEGIN(name)
	#define SMALL_STRAIN_MECH_PROFILE_END()
#endif

using namespace std;

namespace ug {
namespace SmallStrainMechanics{

//////// Volume forces (IMPORT)
template<typename TDomain>
void SmallStrainMechanicsElemDisc<TDomain>::
set_volume_forces(SmartPtr<CplUserData<MathVector<dim>, dim> > user)
{
	m_imVolForce.set_data(user);
}

template<typename TDomain>
void SmallStrainMechanicsElemDisc<TDomain>::
set_volume_forces(number vel_x)
{
	SmartPtr<ConstUserVector<dim> > force(new ConstUserVector<dim>());
	force->set_all_entries(vel_x);
	set_volume_forces(force);
}


template<typename TDomain>
void SmallStrainMechanicsElemDisc<TDomain>::
set_volume_forces(number vel_x, number vel_y)
{
	UG_THROW("SmallStrainMechanicsElemDisc: Setting force vector of dimension 2"
					" to a Discretization for world dim " << dim);
}

template<>
void SmallStrainMechanicsElemDisc<Domain2d>::
set_volume_forces(number vel_x, number vel_y)
{
	SmartPtr<ConstUserVector<dim> > force(new ConstUserVector<dim>());
	force->set_entry(0, vel_x);
	force->set_entry(1, vel_y);
	set_volume_forces(force);
}

template<typename TDomain>
void SmallStrainMechanicsElemDisc<TDomain>::
set_volume_forces(number vel_x, number vel_y, number vel_z)
{
	UG_THROW("SmallStrainMechanicsElemDisc: Setting force vector of dimension 3"
					" to a Discretization for world dim " << dim);
}

template<>
void SmallStrainMechanicsElemDisc<Domain3d>::
set_volume_forces(number vel_x, number vel_y, number vel_z)
{
	SmartPtr<ConstUserVector<dim> > force(new ConstUserVector<dim>());
	force->set_entry(0, vel_x);
	force->set_entry(1, vel_y);
	force->set_entry(2, vel_z);
	set_volume_forces(force);
}


#ifdef UG_FOR_LUA
template<typename TDomain>
void SmallStrainMechanicsElemDisc<TDomain>::
set_volume_forces(const char* fctName)
{
	set_volume_forces(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
}
#endif

/////////////// PRESSURE (IMPORT)

template<typename TDomain>
void SmallStrainMechanicsElemDisc<TDomain>::
set_pressure(SmartPtr<CplUserData<number, dim> > user) 
{
	m_imPressure.set_data(user);
}

template<typename TDomain>
void SmallStrainMechanicsElemDisc<TDomain>::
set_pressure(number val)
{
	set_pressure(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void SmallStrainMechanicsElemDisc<TDomain>::
set_pressure(const char* fctName)
{
	set_pressure(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif
//////////////////////////////////

template<typename TDomain>
void SmallStrainMechanicsElemDisc<TDomain>::update_geo_elem(
		TBaseElem* elem, DimFEGeometry<dim>& geo)
{
	SmartPtr<TDomain> dom = this->domain();

	typedef typename IElemDisc<TDomain>::domain_type::position_accessor_type
			position_accessor_type;
	const position_accessor_type& aaPos = dom->position_accessor();

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get vertices and extract corner coordinates
	const size_t numVertices = elem->num_vertices();
	for (size_t i = 0; i < numVertices; ++i) {
		coCoord[i] = aaPos[elem->vertex(i)];
	};

	//	prepare geometry for type and order
   	try{
		geo.update(elem, &(coCoord[0]), m_lfeID, m_quadOrder);
	}
	UG_CATCH_THROW("SmallStrainMechanics::update_geo_elem:"
					" Cannot update Finite Element Geometry.");
}


template<typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
init_state_variables(const size_t order)
{
	//	call OutputWriter
	//m_spOutWriter->preprocess();
	SMALL_STRAIN_MECH_PROFILE_BEGIN(SmallStrainMechInit_state_variables);

/*#ifdef UG_PARALLEL
	if (pcl::GetProcRank() == 0){
		m_outFile = fopen("sig_eigen.dat", "w");
	}
#else
	m_outFile = fopen("sig_eigen.dat", "w");
#endif*/

	//	TODO: is this necessary here???? Or do I only need num_ip?
	m_order = order;
	m_lfeID = LFEID(LFEID::LAGRANGE, dim, order);
	UG_LOG("\n");
	//	set default quadrature order if not set by user
	if (!m_bQuadOrderUserDef) {
		m_quadOrder = 2 * order + 1;
		UG_LOG("Default QuadOrder is set in 'init_state_variables';"
				"QuadOrder:" << m_quadOrder << "\n");
	}
	UG_LOG("In init_state: Order: " << m_order << " QuadOrder: "
			<< m_quadOrder << "\n");

	//	clears and then attaches the attachments for the internal variables
	//	to the grid
	SmartPtr<TDomain> dom = this->domain();
	typename TDomain::grid_type& grid = *(dom->grid());

	m_spMatLaw->clear_attachments(grid);
	m_spMatLaw->attach_internal_vars(grid);

	//	here the attachment is attached to the whole grid/multigrid
	// 	and it remains attached until the
	//	destructor of SmallStrainMechanicsElemDisc!
	DimFEGeometry<dim> geo;
	typedef typename TDomain::grid_type::template traits<TBaseElem>::iterator
			ElemIter;
	for (ElemIter iter = grid.template begin<TBaseElem> (); iter
			!= grid.template end<TBaseElem> (); iter++)
	{
		TBaseElem* elem = *iter;
		update_geo_elem(elem, geo);

		m_spMatLaw->init_internal_vars(elem, geo.num_ip());

	}//end (ElemIter)

	SMALL_STRAIN_MECH_PROFILE_END();
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
prep_timestep_elem(const number time, const LocalVector& u,
		GeometricObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	call OutputWriter
	if (m_bOutWriter)
		m_spOutWriter->pre_timestep();
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	//	check, if a material law is set
	if (m_spMatLaw.invalid())
		UG_THROW("No material law set in "
				"SmallStrainMechanicsElemDisc::prep_elem_loop \n");

	//	TODO: this only needs to be checked once at the beginning!
	//	-> prepare_setting?!
	if (!m_spMatLaw->is_initialized())
		m_spMatLaw->init();

	//	pass the material law to the output writer
	if ((!m_bMatLawPassedToOutWriter) && m_spOutWriter.valid()){
		m_spOutWriter->material_law(m_spMatLaw);
		m_spOutWriter->quad_order(m_quadOrder);
		m_bMatLawPassedToOutWriter = true;
	}

	//	request geometry
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	//	prepare geometry for type and order
	try{
		geo.update_local(roid, m_lfeID, m_quadOrder);
	}UG_CATCH_THROW("SmallStrainMechanicsElemDisc::prep_elem_loop:"
					" Cannot update Finite Element Geometry.");

	//	set local positions for rhs
	static const int refDim = TElem::dim;
	m_imVolForce.template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), true);
	m_imPressure.template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::fsh_elem_loop()
{}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
prep_elem(const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	request geometry
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	try{
		geo.update(elem, vCornerCoords, m_lfeID, m_quadOrder);
	}
	UG_CATCH_THROW("SmallStrainMechanics::prep_elem:"
					" Cannot update Finite Element Geometry.");

	//	set global positions for rhs
	m_imVolForce.set_global_ips(geo.global_ips(), geo.num_ip());

	//	set pointer to internal variables of elem
	m_spMatLaw->internal_vars(static_cast<TElem*>(elem));
}

//  assemble stiffness jacobian
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u,
		GeometricObject* elem, const MathVector<dim> vCornerCoords[])
{
	SMALL_STRAIN_MECH_PROFILE_BEGIN(SmallStrainMechAddJacA);

	//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	MathMatrix<dim, dim> GradU;
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//	TODO: think about moving the call of 'DisplacementGradient',
		//	since it is not used for Hooke`s linear elastic law!
		m_spMatLaw->template DisplacementGradient<TFEGeom>(GradU, ip, geo, u);
		m_spElastTensor = m_spMatLaw->elasticityTensor(ip, GradU);

		// A) Compute Du:C:Dv = Du:sigma = sigma:Dv
		for (size_t a = 0; a < geo.num_sh(); ++a){ 				// loop shape functions
			for (size_t i = 0; i < (size_t) TDomain::dim; ++i){ // loop component
				for (size_t b = 0; b < geo.num_sh(); ++b){ 		// shape functions
					for (size_t j = 0; j < (size_t) TDomain::dim; ++j) // loop component
					{
						number integrandC = 0.0;

						// Du:C:Dv = Du:sigma = sigma:Dv
						for (size_t K = 0; K < (size_t) dim; ++K){
							for (size_t L = 0; L < (size_t) dim; ++L)
							{
								integrandC += geo.global_grad(ip, a)[K]
										* (*m_spElastTensor)[i][K][j][L]
										* geo.global_grad(ip, b)[L];
							}
						}

						J(i, a, j, b) += integrandC * geo.weight(ip);

						// Du:p*Id = p*Id:Dv goes here..

					} //end (j)
				} //end (b)
			} //end(i)
		} //end(a)
	} //end(ip)
}

//	assemble mass jacobian
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u,
		GeometricObject* elem, const MathVector<dim> vCornerCoords[])
{
	if(m_bAddMassJac)
	{
		//	request geometry
		const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

		for(size_t ip = 0; ip < geo.num_ip(); ++ip){
			for(size_t i = 0; i < geo.num_sh(); ++i){
				for(size_t j = 0; j < geo.num_sh(); ++j)
				{
					// same value for all components
					number value = geo.shape(ip, i) * geo.shape(ip, j) * geo.weight(ip);

					for(size_t c = 0; c < (size_t)dim; ++c){
						J(c, i, c, j) += value;
					}
				}
			}
		}
	}
}

//  assemble stiffness defect
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u,
		GeometricObject* elem, const MathVector<dim> vCornerCoords[])
{
	SMALL_STRAIN_MECH_PROFILE_BEGIN(SmallStrainMechAddDefA);
	//if (0){
		//	request geometry
	/*	static const typename TGeomProvider::Type& geo = TGeomProvider::get();


		for (size_t ip = 0; ip < geo.num_ip(); ++ip)
		{
			MathMatrix<dim, dim> sigma, GradU, eps;

			//	compute cauchy-stress tensor sigma at a ip
			DisplacementGradient<TGeomProvider>(GradU, geo, u, ip);

			for(size_t i = 0; i < (size_t)dim; ++i){
				for(size_t j = 0; j < (size_t)dim; ++j){
					eps[i][j] = 0.5 * (GradU[i][j] + GradU[j][i]);
				}
			}

			// A) Compute Du:C:Dv = Du:sigma = sigma:Dv
			for (size_t a = 0; a < geo.num_sh(); ++a){ // loop shape functions
				for (size_t i = 0; i < (size_t) TDomain::dim; ++i){ // loop component
					for (size_t b = 0; b < geo.num_sh(); ++b){ // shape functions
						for (size_t j = 0; j < (size_t) TDomain::dim; ++j) // loop component
						{
							number integrandC = 0.0;

									// Du:C:Dv = Du:sigma = sigma_Dv
									for (size_t K = 0; K < (size_t) dim; ++K){
										for (size_t L = 0; L < (size_t) dim; ++L)
										{
											integrandC += geo.global_grad(ip, a)[K]
													// * myTensor[i][K][j][L]
													* geo.global_grad(ip, b)[L];
										}
									}

									J(i, a, j, b) += integrandC * geo.weight(ip);
								} //end (j)
							} //end (b)
						} //end(i)
					} //end(a)
		}

		// p* div (v_i)
		for (size_t a = 0; a < geo.num_sh(); ++a)
		{ // loop shapes
			for (size_t i = 0; i < num_fct(); ++i)
			{

				d(i, a) += geo.weight(ip)*0.0;
			}

		}*/
	//}

	//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	MathMatrix<dim, dim> sigma, GradU;
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//	compute cauchy-stress tensor sigma at a ip
		m_spMatLaw->template DisplacementGradient<TFEGeom>(GradU, ip, geo, u);
		m_spMatLaw->stressTensor(sigma, ip, GradU);

		for (size_t a = 0; a < geo.num_sh(); ++a){ // loop shape functions
			for (size_t i = 0; i < num_fct(); ++i) // loop components
			{
				number innerForcesIP = 0.0;

				//	i-th comp. of INTERNAL FORCES at node a:
				for (size_t J = 0; J < (size_t) dim; ++J){
					innerForcesIP += sigma[i][J] * geo.global_grad(ip, a)[J];
				}

				d(i, a) += geo.weight(ip) * innerForcesIP;
			} //end (i)
		} //end (a)
	}//end (ip)

}

//  assemble mass-defect
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u,
		GeometricObject* elem, const MathVector<dim> vCornerCoords[])
{}

//  assemble right-hand-side
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
add_rhs_elem(LocalVector& d, GeometricObject* elem, const MathVector<dim> vCornerCoords[])
{
	// consider volume forces
	if(!m_imVolForce.data_given()) return;

	//	request geometry
	 const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	 for(size_t ip = 0; ip < geo.num_ip(); ++ip){ 	// loop ip
		 for(size_t a = 0; a < geo.num_sh(); ++a){	// loop shape functions
			 for(size_t i = 0; i < num_fct(); ++i){ // loop component
				 d(i,a) += m_imVolForce[ip][i] * geo.shape(ip, a) * geo.weight(ip);
			 }
		 }
	}
}


//	computes the linearized defect w.r.t to the pressure
// $$ \frac{\partial d(u,I_1), ... d(u,I_n)}{\partial I_i} $$
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
lin_def_pressure(const LocalVector& u,
                 std::vector<std::vector<number> > vvvLinDef[],
                 const size_t nip)
{
	//	request geometry
	//static const typename TGeomProvider::Type& geo = TGeomProvider::get();
	static const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
	

	//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//	loop test spaces
		for (size_t i = 0; i < num_fct(); ++i) { // loop component
			for (size_t a = 0; a < geo.num_sh(); ++a) { // loop shape functions


			}
		}
	}
}

//	computes the linearized defect w.r.t to the volume forces
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
lin_def_volume_forces(const LocalVector& u,
	                  std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
	                  const size_t nip)
{
	//	request geometry
		//static const typename TGeomProvider::Type& geo = TGeomProvider::get();
		static const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	//	loop integration points
		for(size_t ip = 0; ip < geo.num_ip(); ++ip)
		{

			for (size_t i = 0; i < num_fct(); ++i)
			{ // loop component

				// 	get ui at integration point
				number shape_ui = 0.0;
				for(size_t a = 0; a < geo.num_sh(); ++a)
					shape_ui += u(i,a) * geo.shape(ip, a);

				for(size_t a = 0; a < geo.num_sh(); ++a)
				 { //	loop test spaces

					//	add to local defect
//					geo.shape(ip, a) * geo.weight(ip);
					VecScale(vvvLinDef[ip][i][a], geo.global_grad(ip, a),
				         	 	 	 	 	 	 geo.weight(ip) * shape_ui);
				}
			}
		}
	}


template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
fsh_timestep_elem(const number time, const LocalVector& u,
		GeometricObject* elem, const MathVector<dim> vCornerCoords[])
{
	SMALL_STRAIN_MECH_PROFILE_BEGIN(SmallStrainMechFshTimeStepElem);

	//	request geometry
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	try{
		geo.update(elem, vCornerCoords, m_lfeID, m_quadOrder);
	}
	UG_CATCH_THROW("SmallStrainMechanics::fsh_timestep_elem:"
					" Cannot update Finite Element Geometry.");

	//  pointer to internal variable of current elem
	m_spMatLaw->internal_vars(static_cast<TElem*>(elem));

	//	call OutputWriter
	if (m_bOutWriter){
		SmartPtr<TDomain> dom = this->domain();
		m_spOutWriter->post_timestep(time, dom, geo, static_cast<TElem*>(elem), u);
	}

	/*MathMatrix<dim, dim> GradU;
	//typedef typename IElemDisc<TDomain>::domain_type::position_accessor_type
	//		position_accessor_type;
	//position_accessor_type& aaPos = dom->position_accessor();
	//	COMPUTE EIGENVALUES OF A STRESSTENSOR AT IP
	if (m_stressEV)
	{
		MathVector<dim> eval_point;
		eval_point[0] = m_eigCoordX;
		if (dim > 1) eval_point[1] = m_eigCoordY;
		if (dim > 2) eval_point[2] = m_eigCoordZ;

		if ((ContainsPoint(static_cast<TElem*>(elem), eval_point, aaPos) == true)
				&& (m_bIP_values_written == false))
		{
			// 	TODO: move next_ips_to_point-call to an instance of
			//	IMechOutputWriter like SmallStrainMechOutput!
			vector<size_t> vNextIP;
			m_spOutWriter->next_ips_to_point(vNextIP, eval_point, geo);

			MathMatrix<dim, dim> SYMMSig, Sig;
			for (vector<size_t>::iterator it = vNextIP.begin();
							it != vNextIP.end(); ++it)
			{
				m_spMatLaw->template DisplacementGradient<TFEGeom>(GradU, *it, geo, u);
				m_spMatLaw->stressTensor(Sig, *it, GradU);

				if (dim == 3)
				{
					MatMultiplyMTM(SYMMSig, Sig);

					MathMatrix<3, 3, number> SigTens;
					for (size_t a = 0; a < 3; ++a)
						for (size_t b = 0; b < 3; ++b){
							SigTens[a][b] = SYMMSig[a][b];
						}

					MathVector<3, number> evMin, evMed, evMax;
					number SiglambdaMin, SiglambdaMed, SiglambdaMax;

					CalculateEigenvalues(SigTens, SiglambdaMin,
							SiglambdaMed, SiglambdaMax, evMin, evMed,
							evMax);

					number sqrSigLambdaMin = sqrt(SiglambdaMin);
					number sqrSigLambdaMed = sqrt(SiglambdaMed);
					number sqrSigLambdaMax = sqrt(SiglambdaMax);

					UG_LOG("minimal EigenValueSigma: " << sqrSigLambdaMin << "\n");
					UG_LOG("medium EigenValueSigma: " << sqrSigLambdaMed << "\n");
					UG_LOG("maximal EigenValueSigma: " << sqrSigLambdaMax << "\n");

					number absSiglambdamax = abs(sqrSigLambdaMax);

				#ifdef UG_PARALLEL
					if (pcl::GetProcRank() == 0){
						fprintf(m_outFile, "%lg %lg \n ", time, absSiglambdamax);
					}
				#else
					fprintf(m_outFile, "%lg %lg \n ", time, absSiglambdamax);
				#endif

				}

				//	to be sure that data is only written once:
				m_bIP_values_written = true;
			} // end vNextIP-iteration
		} //end(if containsPoint)
	} //end(if m_stressEV)*/

	//	COMPUTE NORMAL STRESSES OF THE STRESSTENSOR SIGMA AT IP
	/*if (m_normalStress)
	{
		MathVector<dim> eval_point;
		eval_point[0] = m_normStressCoordX;
		if (dim > 1) eval_point[1] = m_normStressCoordY;
		if (dim > 2) eval_point[2] = m_normStressCoordZ;

		if ((ContainsPoint(static_cast<TElem*>(elem), eval_point, aaPos) == true))
				//&& (m_bIP_values_written == false))
		{
			// 	TODO: move next_ips_to_point-call to an instance of
			//	IMechOutputWriter like SmallStrainMechOutput!
			vector<size_t> vNextIP;
			m_spOutWriter->next_ips_to_point(vNextIP, eval_point, geo);

			MathMatrix<dim, dim> Sig;
			for (vector<size_t>::iterator it = vNextIP.begin();
							it != vNextIP.end(); ++it)
			{
				m_spMatLaw->template DisplacementGradient<TFEGeom>(GradU, *it, geo, u);
				m_spMatLaw->stressTensor(Sig, *it, GradU);

				UG_LOG("At " << geo.global_ip(*it) << ": \n");
				UG_LOG("sigma_xx: " << Sig[0][0] << "\n");
				UG_LOG("sigma_yy: " << Sig[1][1] << "\n");
				UG_LOG("sigma: " << Sig << "\n");
				UG_LOG("######### \n");

				//	to be sure that data is only written once:
				//m_bIP_values_written = true;
			} // end vNextIP-iteration
		} //end(if containsPoint)
	}*/

	//	UPDATE PLASTIC STRAINS 'eps_p' AND HARDENING VARIABLE 'alpha'
	MathMatrix<dim, dim> GradU;
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//  update the internal (plastic) variables strain_p_new and alpha
		m_spMatLaw->template DisplacementGradient<TFEGeom>(GradU, ip, geo, u);
		m_spMatLaw->update_internal_vars(ip, GradU);
	}

}

/*template<typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
normal_stress_strain_loc(LocalVector& locDevSigma, LocalVector& locSigma,
		LocalVector& locEps, TBaseElem* elem, const LocalVector& locU)
{
	SmartPtr<TDomain> dom = this->domain();

//	get all neighbor elems which share a vertex with the given element 'elem'
	typename TDomain::grid_type& grid = *(dom->grid());
	vector<TBaseElem*> vNeighborElems;
	CollectNeighbors(vNeighborElems, elem, grid, NHT_VERTEX_NEIGHBORS);
	typedef typename vector<TBaseElem*>::iterator neighborElemIter;

	typedef typename IElemDisc<TDomain>::domain_type::position_accessor_type
			position_accessor_type;
	position_accessor_type& aaPos = dom->position_accessor();

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get vertices and extract corner coordinates
	const size_t numVertices = elem->num_vertices();
	for (size_t i = 0; i < numVertices; ++i) {
		coCoord[i] = aaPos[elem->vertex(i)];
	};

	//	prepare geometry for type and order
	DimFEGeometry<dim> geo;
   	try{
		geo.update(elem, &(coCoord[0]), m_lfeID, m_quadOrder);
	}
   	UG_CATCH_THROW("SmallStrainMechanicsElemDisc::normal_stress_strain_loc:"
					" Cannot update Finite Element Geometry.");

	//	m_pElemData = &m_aaElemData[static_cast<TBaseElem*>(elem)];
	//  pointer to internal variable of current elem
	m_spMatLaw->internal_vars(elem);

	MathMatrix<dim, dim> GradU, eps, sigma, devSigma, SumEpsIP, SumSigmaIP;
	number normDevSig, SumNormDevSigIP;

	//	init summation-values
	SumEpsIP = 0.0; SumSigmaIP = 0.0; SumNormDevSigIP = 0.0;

	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//	get displacementGradient (GradU) and Cauchy stress-tensor sigma at ip
		//DisplacementGradient<DimFEGeometry<dim> >(GradU, geo, locU, ip);
		m_spMatLaw->template DisplacementGradient<DimFEGeometry<dim> >(GradU, ip, geo, locU);
		m_spMatLaw->stressTensor(sigma, ip, GradU);

		// 	get linearized strain tensor (eps) at ip
		for (size_t i = 0; i < (size_t) dim; ++i)
			for (size_t J = 0; J < (size_t) dim; ++J)
				eps[i][J] = 0.5 * (GradU[i][J] + GradU[J][i]);

		//	get norm of deviator of sigma
		MatDeviatorTrace(sigma, devSigma);
		normDevSig = MatFrobeniusNorm(devSigma);

		//	add values to local vector
		for (size_t i = 0; i < (size_t) dim; ++i)
		{
			SumEpsIP[i][i] += eps[i][i];
			SumSigmaIP[i][i] += sigma[i][i];
		}
		SumNormDevSigIP += normDevSig;

	} //end (ip)

	for (size_t co = 0; co < numVertices; ++co) // loop corner
	{
		//	TODO: is this geometry update necessary here?

		//	update geometry
		try{
			geo.update(elem, &(coCoord[0]), m_lfeID, m_quadOrder);
		}
		UG_CATCH_THROW("SmallStrainMechanicsElemDisc::normal_stress_strain_loc:"
						" Cannot update Finite Element Geometry.");

		//	init 'elemsWithCo' with 1, because co lies in 'elem'
		size_t elemsWithCo = 1;

		//	iterate neighbor elems and count how many elems contain the corner 'co'
		for (neighborElemIter it = vNeighborElems.begin();
				it != vNeighborElems.end(); ++it){
			if (ContainsPoint(*it, coCoord[co], aaPos) == true)
				elemsWithCo++;
		}

		//	scaling factor for averaging values out of all ip�s
		//	of the associated elements of corner co
		size_t scaleFac = elemsWithCo * geo.num_ip();

		//	add values to local vector
		for (size_t i = 0; i < (size_t) dim; ++i) // loop component
		{
			locEps(i, co) += SumEpsIP[i][i] / scaleFac;
			locSigma(i, co) += SumSigmaIP[i][i] / scaleFac;
		}
		locDevSigma(0, co) += SumNormDevSigIP / scaleFac;
	} //end(co)


	//	temporarily output for linear-elastic-benchmark-output:
	for (size_t co = 0; co < numVertices; ++co)
	{
		if ((coCoord[co][0] == 90.0) && (coCoord[co][1] == 0.0))
		{
			UG_LOG("At 2 (90.0,0.0): \n");
			for (size_t i = 0; i < (size_t) dim; ++i) {
				UG_LOG("u("<< i << "," << co << "): " << locU(i,co) << "\n");
			}
			UG_LOG("\n");

			vector<size_t> vNextIP;
			m_spOutWriter->next_ips_to_point(vNextIP, coCoord[co], geo);

			MathMatrix<dim, dim> sigma;
			for (vector<size_t>::iterator it = vNextIP.begin();
							it != vNextIP.end(); ++it)
			{
				//	get displacementGradient (GradU) at ip
				m_spMatLaw->template DisplacementGradient<DimFEGeometry<dim> >(GradU, *it, geo, locU);
				//	get the Cauchy Stress Tensor (sigma) at ip
				m_spMatLaw->stressTensor(sigma, *it, GradU);

				UG_LOG("At " << geo.global_ip(*it) << ": \n");
				UG_LOG("sigma_yy: " << sigma[1][1] << "\n");
				UG_LOG("sigma: " << sigma << "\n");
			}
		}

		if ((coCoord[co][0] == 100.0) && (coCoord[co][1] == 100.0))
		{
			UG_LOG("At 4 (100.0,100.0): \n");
			for (size_t i = 0; i < (size_t) dim; ++i) {
				UG_LOG("u("<< i << "," << co << "): " << locU(i,co) << "\n");
			}
			UG_LOG("\n");
		}

		if ((coCoord[co][0] == 0.0) && (coCoord[co][1] == 100.0))
		{
			UG_LOG("At 5 (0.0,100.0): \n");
			for (size_t i = 0; i < (size_t) dim; ++i) {
				UG_LOG("u("<< i << "," << co << "): " << locU(i,co) << "\n");
			}
			UG_LOG("\n");
		}

	};
}*/


/*template <typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
close_gnuplot_file()
{}*/

template <typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
set_assemble_funcs()
{
	//	set default quadrature order if not set by user
	if (!m_bQuadOrderUserDef) {
		m_quadOrder = 2 * m_order + 1;
	}
	//	set all non-set orders
	else
	{
		if (m_quadOrder < 0){
			m_quadOrder = 2 * m_order + 1;
		}
	}

	register_all_fe_funcs(m_order, m_quadOrder);

}

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
SmallStrainMechanicsElemDisc<TDomain>::
SmallStrainMechanicsElemDisc(const char* functions, const char* subsets) :
			IElemDisc<TDomain> (functions, subsets),
			m_spMatLaw(SPNULL), m_spElastTensor(SPNULL), m_spOutWriter(SPNULL),
			m_bMatLawPassedToOutWriter(false), m_bAddMassJac(false),
			m_bOutWriter(false)
			//,m_stressEV(false), m_normalStress(false),
			//m_bIP_values_written(false)
{
	//	check number of functions
	if (this->num_fct() != (size_t) dim)
		UG_THROW("Wrong number of functions: The ElemDisc 'SmallStrainMechanics'"
				" needs exactly "<<dim<<" symbolic function.");

	//	register imports
	this->register_import(m_imPressure);
	this->register_import(m_imVolForce);
	m_imVolForce.set_rhs_part();

	//	set defaults
	m_order = 1;
	m_bQuadOrderUserDef = false;
	m_quadOrder = -1;

	//	update assemble functions
	set_assemble_funcs();

	// TODO: check, if OutWriter is set (by means of a bool var)
	// if not set m_spOutWriter to a default instance of IMechOutputWriter!
}


template <typename TDomain>
SmallStrainMechanicsElemDisc<TDomain>::
~SmallStrainMechanicsElemDisc()
{
	//if(m_bOutWriter)
	//	m_spOutWriter->postprocess();

	SmartPtr<TDomain> dom = this->domain();
	typename TDomain::grid_type& grid = *(dom->grid());
	m_spMatLaw->clear_attachments(grid);
}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void SmallStrainMechanicsElemDisc<Domain1d>::register_all_fe_funcs(int order,
		int quadOrder)
{
	//	Edge
	register_fe_func<Edge, DimFEGeometry<dim> > ();
}
#endif

#ifdef UG_DIM_2
template<>
void SmallStrainMechanicsElemDisc<Domain2d>::register_all_fe_funcs(int order,
		int quadOrder)
{
	if (quadOrder != 2 * order + 1) {
		register_fe_func<Triangle, DimFEGeometry<dim> > ();
		register_fe_func<Quadrilateral, DimFEGeometry<dim> > ();
	}

	//	special compiled cases

	//	Triangle
	switch (order)
	{
		case 1:{
			typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 1> ,
					GaussQuadrature<ReferenceTriangle, 3> > FEGeom;
			register_fe_func<Triangle, FEGeom > (); break;
		}
		case 2:{
			typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 2> ,
					GaussQuadrature<ReferenceTriangle, 5> > FEGeom;
			register_fe_func<Triangle, FEGeom > (); break;
		}
		case 3:{
			typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 3> ,
					GaussQuadrature<ReferenceTriangle, 7> > FEGeom;
			register_fe_func<Triangle, FEGeom > (); break;
		}
		default: register_fe_func<Triangle, DimFEGeometry<dim> > (); break;
	}

	//	Quadrilateral
	switch (order)
	{
		case 1:{
			typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<
					ReferenceQuadrilateral, 1> , GaussQuadrature<
					ReferenceQuadrilateral, 3> > FEGeom;
			register_fe_func<Quadrilateral, FEGeom > (); break;
		}
		case 2:{
			typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<
					ReferenceQuadrilateral, 2> , GaussQuadrature<
					ReferenceQuadrilateral, 7> > FEGeom;
			register_fe_func<Quadrilateral, FEGeom > (); break;
		}
		case 3:{
			typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<
					ReferenceQuadrilateral, 3> , GaussQuadrature<
					ReferenceQuadrilateral, 11> > FEGeom;
			register_fe_func<Quadrilateral, FEGeom > (); break;
		}
		default: register_fe_func<Quadrilateral, DimFEGeometry<dim> > (); break;
	}
}
#endif

#ifdef UG_DIM_3
template<>
void SmallStrainMechanicsElemDisc<Domain3d>::register_all_fe_funcs(int order,
		int quadOrder)
{
	if (quadOrder != 2 * order + 1) {
		register_fe_func<Tetrahedron, DimFEGeometry<dim> > ();
		register_fe_func<Prism, DimFEGeometry<dim> > ();
		register_fe_func<Pyramid, DimFEGeometry<dim> > ();
		register_fe_func<Hexahedron, DimFEGeometry<dim> > ();
	}

	//	special compiled cases

	//	Tetrahedron
	switch (order) 
	{
		case 1:{
			typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<
				ReferenceTetrahedron, 1> , GaussQuadrature<
				ReferenceTetrahedron, 3> > FEGeom;
			register_fe_func<Tetrahedron, FEGeom > (); break;
		}
		case 2:{
			typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<
					ReferenceTetrahedron, 2> , GaussQuadrature<
					ReferenceTetrahedron, 5> > FEGeom;
			register_fe_func<Tetrahedron, FEGeom > (); break;
		}
		case 3:{
			typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<
					ReferenceTetrahedron, 3> , GaussQuadrature<
					ReferenceTetrahedron, 7> > FEGeom;
			register_fe_func<Tetrahedron, FEGeom > (); break;
		}
		default: register_fe_func<Tetrahedron, DimFEGeometry<dim> > (); break;

	}

	//	Prism
	switch (order)
	{
		case 1:{
			typedef FEGeometry<Prism, dim, LagrangeLSFS<ReferencePrism, 1> ,
					GaussQuadrature<ReferencePrism, 2> > FEGeom;
			register_fe_func<Prism, FEGeom > (); break;
		}
		default:
			register_fe_func<Prism, DimFEGeometry<dim> > (); break;
	}

	//	Pyramid
	switch (order)
	{
		default: register_fe_func<Pyramid, DimFEGeometry<dim> > (); break;
	}

	//	Hexahedron
	switch (order)
	{
		case 1:{
			if (quadOrder == 2){
				typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<
				ReferenceHexahedron, 1> , GaussQuadrature<
				ReferenceHexahedron, 2> > FEGeom;
				register_fe_func<Hexahedron, FEGeom > (); break;
			}
			else{
				typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<
				ReferenceHexahedron, 1> , GaussQuadrature<
				ReferenceHexahedron, 3> > FEGeom;
				register_fe_func<Hexahedron, FEGeom > (); break;
			}
		}
		case 2:{
			typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<
					ReferenceHexahedron, 2> , GaussQuadrature<
					ReferenceHexahedron, 7> > FEGeom;
			register_fe_func<Hexahedron, FEGeom > (); break;
		}
		case 3:{
			typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<
					ReferenceHexahedron, 3> , GaussQuadrature<
					ReferenceHexahedron, 11> > FEGeom;
			register_fe_func<Hexahedron, FEGeom > (); break;
		}
		default: register_fe_func<Hexahedron, DimFEGeometry<dim> > (); break;
	}
}
#endif

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::register_fe_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_add_elem(true);
	this->set_prep_timestep_elem_fct(id,
			&T::template prep_timestep_elem<TElem, TFEGeom>);
	this->set_prep_elem_loop_fct(id,
			&T::template prep_elem_loop<TElem, TFEGeom>);
	this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFEGeom>);
	this->set_fsh_elem_loop_fct(id,
					&T::template fsh_elem_loop<TElem, TFEGeom>);

	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFEGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFEGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFEGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFEGeom>);
	this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFEGeom>);

	this->set_fsh_timestep_elem_fct(id,
			&T::template fsh_timestep_elem<TElem, TFEGeom>);

	//	set computation of linearized defect w.r.t velocity
	m_imVolForce.set_fct(id, this, &T::template lin_def_volume_forces<TElem, TFEGeom>);
	m_imPressure.set_fct(id, this, &T::template lin_def_pressure<TElem, TFEGeom>);

	//	exports
//	m_exValue->	  template set_fct<T,refDim>(id, this, &T::template ex_value_fe<TElem, TGeomProvider>);
//	m_exGrad->    template set_fct<T,refDim>(id, this, &T::template ex_grad_fe<TElem, TGeomProvider>);

}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////


#ifdef UG_DIM_2
template class SmallStrainMechanicsElemDisc<Domain2d> ;
#endif
#ifdef UG_DIM_3
template class SmallStrainMechanicsElemDisc<Domain3d> ;
#endif

} // namespace SmallStrainMechanics
} // namespace ug

