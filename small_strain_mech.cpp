/*
 * small_strain_mech.cpp
 *
 *  Created on: 16.05.2012
 *      Author: raphaelprohl, Andreas Vogel
 */

#include "small_strain_mech.h"

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

//////// Volume forces (IMPORT)
template<typename TDomain>
void SmallStrainMechanicsElemDisc<TDomain>::
set_viscous_forces(SmartPtr<CplUserData<MathVector<dim>, dim> > user1, SmartPtr<CplUserData<MathVector<dim>, dim> > user2)
{
	m_imViscousForces[0].set_data(user1);
	m_imViscousForces[1].set_data(user2);
}

//////////////////////////////////

template<typename TDomain>
void SmallStrainMechanicsElemDisc<TDomain>::
update_geo_elem(TBaseElem* elem, DimFEGeometry<dim>& geo)
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
	SMALL_STRAIN_MECH_PROFILE_BEGIN(SmallStrainMechInit_state_variables);

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
		GridObject* elem, const MathVector<dim> vCornerCoords[])
{
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
	if (m_bOutWriter && (!m_bMatLawPassedToOutWriter)){
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
	m_imViscousForces[0].template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), true);
	m_imViscousForces[1].template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), true);
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::fsh_elem_loop()
{}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
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
	m_imViscousForces[0].set_global_ips(geo.global_ips(), geo.num_ip());
	m_imViscousForces[1].set_global_ips(geo.global_ips(), geo.num_ip());
	m_imPressure.set_global_ips(geo.global_ips(), geo.num_ip());


	//	set pointer to internal variables of elem
	m_spMatLaw->internal_vars(static_cast<TElem*>(elem));
}

//  assemble stiffness jacobian
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u,
		GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	SMALL_STRAIN_MECH_PROFILE_BEGIN(SmallStrainMechAddJacA);

	//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	if (m_spMatLaw->elastTensIsConstant())
	{
		//	material law with constant elasticity tensors
		m_spElastTensor = m_spMatLaw->elasticityTensor();
	}

	MathMatrix<dim, dim> GradU;
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		if (!(m_spMatLaw->elastTensIsConstant()))
		{
			//	material law with non constant elasticity tensors
			m_spMatLaw->template DisplacementGradient<TFEGeom>(GradU, ip, geo, u);
			m_spElastTensor = m_spMatLaw->elasticityTensor(ip, GradU);
		}

		// A) Compute Du:C:Dv = Du:sigma = sigma:Dv
		for (size_t a = 0; a < geo.num_sh(); ++a) 				// loop shape functions
			for (size_t i = 0; i < (size_t) TDomain::dim; ++i) // loop component
				for (size_t b = 0; b < geo.num_sh(); ++b) 		// shape functions
					for (size_t j = 0; j < (size_t) TDomain::dim; ++j) // loop component
					{
						number integrandC = 0.0;

						// Du:C:Dv = Du:sigma = sigma:Dv
						for (size_t K = 0; K < (size_t) dim; ++K)
							for (size_t L = 0; L < (size_t) dim; ++L)
							{
								integrandC += geo.global_grad(ip, a)[K]
										* (*m_spElastTensor)[i][K][j][L]
										* geo.global_grad(ip, b)[L];
							}

						J(i, a, j, b) += integrandC * geo.weight(ip);

						// Du:p*Id = p*Id:Dv goes here..

					} //end (j)

	} //end(ip)
}

//	assemble mass jacobian
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u,
		GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	if(m_spMatLaw->needs_to_add_jac_m())
	{
		//	request geometry
		const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

		for(size_t ip = 0; ip < geo.num_ip(); ++ip)
			for(size_t i = 0; i < geo.num_sh(); ++i)
				for(size_t j = 0; j < geo.num_sh(); ++j)
				{
					// same value for all components
					number value = geo.shape(ip, i) * geo.shape(ip, j) * geo.weight(ip);

					for(size_t c = 0; c < (size_t)dim; ++c)

					{
						J(c, i, c, j) += value*m_massScale;

					}
				}
	}
}

//  assemble stiffness defect
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u,
		GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	SMALL_STRAIN_MECH_PROFILE_BEGIN(SmallStrainMechAddDefA);

	//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	// a) Cauchy tensor
	MathMatrix<dim, dim> sigma, GradU;
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//	compute cauchy-stress tensor sigma at a ip
		m_spMatLaw->template DisplacementGradient<TFEGeom>(GradU, ip, geo, u);
		m_spMatLaw->stressTensor(sigma, ip, GradU);

		for (size_t sh = 0; sh < geo.num_sh(); ++sh) // loop shape functions
			for (size_t i = 0; i < num_fct(); ++i) // loop components
			{
				number innerForcesIP = 0.0;

				//	i-th comp. of INTERNAL FORCES at node a:
				for (size_t j = 0; j < (size_t) dim; ++j)
					innerForcesIP += sigma[i][j] * geo.global_grad(ip, sh)[j];

				d(i, sh) += geo.weight(ip) * innerForcesIP;
			} //end (i)

	}//end (ip)




}

//  assemble mass-defect
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u,
		GridObject* elem, const MathVector<dim> vCornerCoords[])
{}

//  assemble right-hand-side d(i,sh)
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	// a) volume forces: $$  \vec F*Phi =(grad p)*Phi$$
	if(m_imVolForce.data_given()) {
		// loop ip
		for(size_t ip = 0; ip < geo.num_ip(); ++ip)
		{
			// loop shape
			for(size_t sh = 0; sh < geo.num_sh(); ++sh)
			{
				// loop component
				for(size_t i = 0; i < num_fct(); ++i)
				{
					d(i,sh) += geo.weight(ip) * geo.shape(ip, sh) * m_imVolForce[ip][i];
				}
			}
		}
	}

	// b) scalar contribution, e.g, p * sum dx_i Phi_sh,i
	if(m_imPressure.data_given()) {

		for(size_t sh = 0; sh < geo.num_sh(); ++sh)
		{	// shape functions (sh=ux, uy, uz)
			for(size_t i = 0; i < num_fct(); ++i)
			{ //  component (i=1,2,3)
				for(size_t ip = 0; ip < geo.num_ip(); ++ip)
				{ 	// loop ip
					number divUip = geo.global_grad(ip, sh)[i];
					d(i,sh) += geo.weight(ip)*divUip*m_imPressure[ip];
				}
			}
		}
	} // pressure data



/*
	// b) pressure part: p(ip) div (V_sh,ip)
	// if (m_imPressure.data_given()) {
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{ // loop IPs

		for (size_t sh = 0; sh < geo.num_sh(); ++sh)
		{ // loop shape functions

			for (size_t i = 0; i < num_fct(); ++i)
			{ // loop components

				number divIP = 0.0;
				for (size_t j = 0; j < (size_t) dim; ++j)
				{
					divIP += geo.global_grad(ip, sh)[j];
				}
				d(i, sh) += geo.weight(ip)*divIP; //import(ip)

			} // end(i)
		}// end (sh)
	} // end (ip)
*/

	// c) viscous stresses:
	// $$ v0^T (grad Phi + grad Phi^T) v1 = v0^T (grad Phi) v1 + ((grad Phi) v0)^T v1 $$

	if (m_imViscousForces[0].data_given() &&
			m_imViscousForces[1].data_given())
	{
		for (size_t ip = 0; ip < geo.num_ip(); ++ip)
		{ // loop IPs

			for (size_t sh = 0; sh < geo.num_sh(); ++sh)
			{ // loop shape functions

				number gradPhi_v0 = VecDot(geo.global_grad(ip, sh), m_imViscousForces[0][ip]);
				number gradPhi_v1 = VecDot(geo.global_grad(ip, sh), m_imViscousForces[1][ip]);
				/*for (size_t i = 0; i < num_fct(); ++i)
				{
					// (grad Phi) v0
					gradPhi_v0 += geo.global_grad(ip, sh)[i]*m_imViscousForces[0][ip][i];
					// (grad Phi) v1
					gradPhi_v1 += geo.global_grad(ip, sh)[i]*m_imViscousForces[1][ip][i];
				}
*/
				for (size_t i = 0; i < num_fct(); ++i)
				{
					d(i, sh) += geo.weight(ip)*m_imViscousForces[0][ip][i]*gradPhi_v1;
					d(i, sh) += geo.weight(ip)*m_imViscousForces[1][ip][i]*gradPhi_v0;

				} // end(i)
			}// end (sh)
		} // end (ip)
	}

}


/**	computes the linearized defect w.r.t to the pressure:
	$$ \frac{\partial d(u,I_1), ... d(u,I_n)}{\partial I_i} $$
	(since d is a vector, we obtain a bunch of vectors)
 */
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
	UG_ASSERT(nip == geo.num_ip(), "Number of integration points does not match!");

	for(size_t sh = 0; sh < geo.num_sh(); ++sh)
	{ // loop test spaces
		for (size_t i = 0; i < num_fct(); ++i)
		{ // loop component
			//	loop integration points
			for(size_t ip = 0; ip < geo.num_ip(); ++ip)
			{
				vvvLinDef[ip][i][sh] = geo.weight(ip)*geo.global_grad(ip, sh)[i];
			}
		}
	}

}

/**	computes the linearized defect w.r.t to the volume forces
	$$ \frac{\partial d(u,I_1), ... d(u,I_n)}{\partial I_i} $$
	(since d is a vector, we obtain a bunch of diagonal matrices)
*/
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
lin_def_volume_forces(const LocalVector& u,
	                  std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
	                  const size_t nip)
{
	//	request geometry
	static const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
	UG_ASSERT(nip == geo.num_ip(), "Number of integration points does not match!");

	//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		for (size_t i = 0; i < num_fct(); ++i)
		{ // loop component

			for(size_t sh = 0; sh < geo.num_sh(); ++sh)
			{ // loop test spaces

				//	add to local defect
				vvvLinDef[ip][i][sh] = 0.0;
				(vvvLinDef[ip][i][sh])[i] = geo.weight(ip)*geo.shape(ip, sh);
				/*VecScale(vvvLinDef[ip][i][sh], geo.global_grad(ip, sh),
						geo.weight(ip) * shape_ui);*/
			}
		}
	}
}

/**	computes the linearized defect w.r.t to the volume forces
	$$ \frac{\partial d(u,I_1), ... d(u,I_n)}{\partial I_i} $$
	(since d is a vector, we obtain a bunch of diagonal matrices)
*/
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
lin_def_viscous_forces0(const LocalVector& u,
	                  std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
	                  const size_t nip)
{
	//	request geometry
	static const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
	UG_ASSERT(nip == geo.num_ip(), "Number of integration points does not match!");

	// loop IPs
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		// loop shape functions
		for (size_t sh = 0; sh < geo.num_sh(); ++sh)
		{
			//number gradPhi_v0 = VecDot(geo.global_grad(ip, sh), m_imViscousForces[0][ip]);
			number gradPhi_v1 = VecDot(geo.global_grad(ip, sh), m_imViscousForces[1][ip]);

			// loop components
			for (size_t i = 0; i < num_fct(); ++i)
			{
				VecScale(vvvLinDef[ip][i][sh], m_imViscousForces[1][ip], geo.weight(ip)*geo.global_grad(ip, sh)[i]);
				(vvvLinDef[ip][i][sh])[i] += geo.weight(ip)*gradPhi_v1;
			} // end(i)
		}// end (sh)
	} // end (ip)end (ip)
}

/**	computes the linearized defect w.r.t to the volume forces
	$$ \frac{\partial d(u,I_1), ... d(u,I_n)}{\partial I_i} $$
	(since d is a vector, we obtain a bunch of diagonal matrices)
*/
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
lin_def_viscous_forces1(const LocalVector& u,
	                  std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
	                  const size_t nip)
{
	//	request geometry
		static const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
		UG_ASSERT(nip == geo.num_ip(), "Number of integration points does not match!");

		// loop IPs
		for (size_t ip = 0; ip < geo.num_ip(); ++ip)
		{
			// loop shape functions
			for (size_t sh = 0; sh < geo.num_sh(); ++sh)
			{
				number gradPhi_v0 = VecDot(geo.global_grad(ip, sh), m_imViscousForces[0][ip]);
				//number gradPhi_v1 = VecDot(geo.global_grad(ip, sh), m_imViscousForces[1][ip]);

				// loop components
				for (size_t i = 0; i < num_fct(); ++i)
				{
					VecScale(vvvLinDef[ip][i][sh], m_imViscousForces[0][ip], geo.global_grad(ip, sh)[i]);
					(vvvLinDef[ip][i][sh])[i] += gradPhi_v0;
				} // end(i)
			}// end (sh)
		} // end (ip)end (ip)

}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
fsh_timestep_elem(const number time, const LocalVector& u,
		GridObject* elem, const MathVector<dim> vCornerCoords[])
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
		m_spOutWriter->post_timestep_elem(time, dom, geo, static_cast<TElem*>(elem), u);
	}

	//	UPDATE PLASTIC STRAINS 'eps_p' AND HARDENING VARIABLE 'alpha'
	MathMatrix<dim, dim> GradU;
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//  update the internal (plastic) variables strain_p_new and alpha
		m_spMatLaw->template DisplacementGradient<TFEGeom>(GradU, ip, geo, u);
		m_spMatLaw->update_internal_vars(ip, GradU);
	}
}


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
			m_bOutWriter(false), m_bMatLawPassedToOutWriter(false),
			m_exDisplacement(new DataExport<MathVector<dim>, dim>(functions)),
			//m_exDisplacement(new VectorDataExport<dim>(functions)),
			m_exDivergence(new DataExport<number, dim>(functions))
{
	//	check number of functions
	if (this->num_fct() != (size_t) dim)
		UG_THROW("Wrong number of functions: The ElemDisc 'SmallStrainMechanics'"
				" needs exactly "<<dim<<" symbolic function.");

	//	register imports
	this->register_import(m_imPressure);
	this->register_import(m_imVolForce);
	this->register_import(m_imViscousForces[0]);
	this->register_import(m_imViscousForces[1]);

	m_imPressure.set_rhs_part();
	m_imVolForce.set_rhs_part();
	m_imViscousForces[0].set_rhs_part();
	m_imViscousForces[1].set_rhs_part();

	//	set defaults
	m_order = 1;
	m_bQuadOrderUserDef = false;
	m_quadOrder = -1;
	m_massScale = 1.0;

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
//	 export functions
////////////////////////////////////////////////////////////////////////////////

///	computes the displacement
template <typename TDomain>
template <typename TElem, typename TFEGeom>
void  SmallStrainMechanicsElemDisc<TDomain>::
ex_displacement_fe(MathVector<dim> vValue[],
			const MathVector<dim> vGlobIP[],
			number time, int si,
			const LocalVector& u,
			GridObject* elem,
			const MathVector<dim> vCornerCoords[],
			const MathVector<TFEGeom::dim> vLocIP[],
			const size_t nip,
			bool bDeriv,
			std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
	//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

	//	reference dimension
	static const int refDim = reference_element_traits<TElem>::dim;

	//	reference object id
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

	//	FE
	if(vLocIP == geo.local_ips())
	{
		//	Loop ip
		for(size_t ip = 0; ip < geo.num_ip(); ++ip)
		{
			MathVector<dim> &uip=vValue[ip];
			VecSet(uip, 0.0);

			// loop shape-fcts & components
			for(size_t sh = 0; sh < geo.num_sh(); ++sh)
				for (size_t alpha = 0; alpha < (size_t) dim; ++alpha)
					uip[alpha] += geo.shape(ip, sh) * u(alpha,sh);

			if(bDeriv)
			{
				for(size_t sh = 0; sh < geo.num_sh(); ++sh)
					for (size_t alpha = 0; alpha < (size_t) dim; ++alpha)
					{
						vvvDeriv[ip][alpha][sh] = 0.0;
						vvvDeriv[ip][alpha][sh][alpha] = geo.shape(ip, sh);
					}
			}
		}
	}
	// 	general case
	else
	{
		// 	general case
		try{
			//	request for trial space
			const LocalShapeFunctionSet<refDim>& rTrialSpace
			= LocalFiniteElementProvider::get<refDim>(roid, m_lfeID);

			//	number of shape functions
			const size_t numSH = rTrialSpace.num_sh();

			//	storage for values of shape functions at ip
			std::vector<number> vShape(numSH);

			//	loop ips
			for(size_t ip = 0; ip < nip; ++ip)
			{
				//	evaluate at shapes at ip
				rTrialSpace.shapes(vShape, vLocIP[ip]);

				MathVector<dim> &uip=vValue[ip];
				VecSet(uip, 0.0);

				// loop shape-fcts & components
				for(size_t sh = 0; sh < numSH; ++sh)
					for (size_t alpha = 0; alpha < (size_t) dim; ++alpha){
						double uash = u(alpha,sh);
						uip[alpha] += uash * vShape[sh];
					}


				if(bDeriv)
				{
					for(size_t sh = 0; sh < numSH; ++sh)
						for (size_t alpha = 0; alpha < (size_t) dim; ++alpha)
						{ vvvDeriv[ip][alpha][sh] = 0.0;
							vvvDeriv[ip][alpha][sh][alpha] = vShape[sh]; }
				}

			}

		}
		UG_CATCH_THROW("SmallStrainMechanicsElemDisc::ex_displacement_fe: trial space missing.");
	}
}

///	computes the divergence of displacement
template <typename TDomain>
template <typename TElem, typename TFEGeom>
void SmallStrainMechanicsElemDisc<TDomain>::
ex_divergence_fe(number vValue[],
						const MathVector<dim> vGlobIP[],
						number time, int si,
						const LocalVector& u,
						GridObject* elem,
						const MathVector<dim> vCornerCoords[],
						const MathVector<TFEGeom::dim> vLocIP[],
						const size_t nip,
						bool bDeriv,
						std::vector<std::vector<number> > vvvDeriv[])
{
	//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

	//	reference dimension
	static const int refDim = reference_element_traits<TElem>::dim;

	//	reference object id
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

	//	FE
	if(vLocIP == geo.local_ips())
	{
		//	Loop ip
		for(size_t ip = 0; ip < geo.num_ip(); ++ip)
		{
			number &divip=vValue[ip];
			divip = 0.0;

			// loop shape-fcts & components
			for(size_t sh = 0; sh < geo.num_sh(); ++sh)
				for (size_t alpha = 0; alpha < (size_t) dim; ++alpha)
				{ divip += geo.global_grad(ip, sh)[alpha] * u(alpha,sh); }

			if(bDeriv){
				for(size_t sh = 0; sh < geo.num_sh(); ++sh)
					for (size_t alpha = 0; alpha < (size_t) dim; ++alpha)
					{
						vvvDeriv[ip][alpha][sh] = geo.global_grad(ip, sh)[alpha];
					}
			}
		}
	}

	else
	{
		// 	general case
		try{
			//	request for trial space
			const LocalShapeFunctionSet<refDim>& rTrialSpace
			= LocalFiniteElementProvider::get<refDim>(roid, m_lfeID);

			//	number of shape functions
			const size_t numSH = rTrialSpace.num_sh();

			//	storage for gradient function at ip
			std::vector<MathVector<refDim> > vLocGrad(numSH);
			MathVector<refDim> locGradA, globGradA;

			//	Reference Mapping
			MathMatrix<dim, refDim> JTInv;
			ReferenceMapping<ref_elem_type, dim> mapping(vCornerCoords);

			//	loop ips
			for(size_t ip = 0; ip < nip; ++ip)
			{
				// computing divergence:
				number &divip=vValue[ip];
				divip = 0.0;

				//	evaluate shape gradients at ip
				rTrialSpace.grads(vLocGrad, vLocIP[ip]);
				mapping.jacobian_transposed_inverse(JTInv, vLocIP[ip]);


				// loop components
				for (size_t a = 0; a < (size_t) dim; ++a)
				{
					//	compute grad at ip
					VecSet(locGradA, 0.0);
					for(size_t sh = 0; sh < numSH; ++sh)
					{
						VecScaleAppend(locGradA, u(a, sh), vLocGrad[sh]);
					}

					//	compute global grad
					MatVecMult(globGradA, JTInv, locGradA);

					// contribution to divergence
					divip += globGradA[a];

				}

				if(bDeriv){
					assert(0);
					for(size_t sh = 0; sh < geo.num_sh(); ++sh)
						for (size_t alpha = 0; alpha < (size_t) dim; ++alpha)
						{
							vvvDeriv[ip][alpha][sh] /*= 0.0;
							vvvDeriv[ip][alpha][sh][alpha] */= geo.global_grad(ip, sh)[alpha];
						}
				}
			}

		}
		UG_CATCH_THROW("SmallStrainMechanicsElemDisc::ex_divergence_fe: trial space missing.");

	}
}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void SmallStrainMechanicsElemDisc<Domain1d>::register_all_fe_funcs(int order,
		int quadOrder)
{
	//	RegularEdge
	register_fe_func<RegularEdge, DimFEGeometry<dim> > ();
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
	static const int refDim = reference_element_traits<TElem>::dim;

	this->clear_add_fct(id);

	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFEGeom>);
	this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFEGeom>);
	this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFEGeom>);

	this->set_prep_timestep_elem_fct(id, &T::template prep_timestep_elem<TElem, TFEGeom>);
	this->set_fsh_timestep_elem_fct(id, &T::template fsh_timestep_elem<TElem, TFEGeom>);

	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFEGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFEGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFEGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFEGeom>);
	this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFEGeom>);

	// specify computation of linearized defect w.r.t imports
	m_imVolForce.set_fct(id, this, &T::template lin_def_volume_forces<TElem, TFEGeom>);
	m_imPressure.set_fct(id, this, &T::template lin_def_pressure<TElem, TFEGeom>);
	m_imViscousForces[0].set_fct(id, this, &T::template lin_def_viscous_forces0<TElem, TFEGeom>);
	m_imViscousForces[1].set_fct(id, this, &T::template lin_def_viscous_forces1<TElem, TFEGeom>);


	// exports
	m_exDisplacement->  template set_fct<T,refDim>(id, this, &T::template ex_displacement_fe<TElem, TFEGeom>);
	m_exDivergence->    template set_fct<T,refDim>(id, this, &T::template ex_divergence_fe<TElem, TFEGeom>);
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

