/*
 * small_strain_mech_plugin.cpp
 *
 *  Created on: 16.05.2012
 *      Author: raphaelprohl, Andreas Vogel
 */

#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "lib_disc/function_spaces/grid_function.h"

#include "small_strain_mech.h"
#include "small_strain_mech_output.h"
#include "contact/contact.h"

#include "material_laws/hooke_law.h"
#include "material_laws/mat_law_interface.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace SmallStrainMechanics{

/** 
 *  \defgroup small_strain_mechanics Small Strain Mechanics
 *  \ingroup plugins_core
 *  \{
 */

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts.
 * All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

//	typedef
	typedef GridFunction<TDomain, TAlgebra> function_type;

//	SmallStrainMechanics Output
	{
		typedef SmallStrainMechOutput<TDomain, function_type> T;
		string name = string("SmallStrainMechOutput").append(suffix);
		reg.add_class_<T>(name, grp)
			.template add_constructor<void (*)(SmartPtr<SmallStrainMechanicsElemDisc<TDomain> >)>("domain disc")
			.add_method("normal_stresses_strains", &T::normal_stresses_strains, "", "sigma#epsilon#u",
					"computes stress tensor sigma, the deviatoric part of sigma and linearized strain tensor epsilon")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SmallStrainMechOutput", tag);
	}

//	Contact Disc for SmallStrainMechanics-contact problems
	{
		typedef ContactSmallStrainMechanics<TDomain, function_type> T;
		typedef ILagrangeMultiplierDisc<TDomain, function_type> TBase;
		string name = string("ContactSmallStrainMechanics").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<SmallStrainMechanicsElemDisc<TDomain> >)>("domain disc")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ContactSmallStrainMechanics", tag);
	}

}

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	static const int dim = TDomain::dim;
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

//	SmallStrainMechanics (i.e. problems of 'Linear Elasticity'
//	or of 'Linear Elasticity + a plasticity for infinitesimal strains')
	{
		typedef SmallStrainMechanicsElemDisc<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("SmallStrainMechanics").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function#Subsets")
			.add_method("set_elasticity_tensor_orthotropic", &T::set_elasticity_tensor_orthotropic, "", "C11#C12#C13#C22#C23#C33#C44#C55#C66")
			.add_method("set_hooke_elasticity_tensor", &T::set_hooke_elasticity_tensor, "", "lambda#mu")
			.add_method("set_hooke_elasticity_tensor_E_nu", &T::set_hooke_elasticity_tensor_E_nu, "", "E#nu")
			.add_method("set_quad_order", &T::set_quad_order, "", "order")

			.add_method("set_volume_forces", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_volume_forces),"", "Force field")
			.add_method("set_volume_forces", static_cast<void (T::*)(number)>(&T::set_volume_forces), "", "F")
			.add_method("set_volume_forces", static_cast<void (T::*)(number,number)>(&T::set_volume_forces), "", "F_x, F_y")
			.add_method("set_volume_forces", static_cast<void (T::*)(number,number,number)>(&T::set_volume_forces), "", "F_x, F_y, F_z")
#ifdef UG_FOR_LUA
			.add_method("set_volume_forces", static_cast<void (T::*)(const char*)>(&T::set_volume_forces), "", "Force field")
#endif

			.add_method("set_pressure", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_pressure), "", "Pressure")
			.add_method("set_pressure", static_cast<void (T::*)(number)>(&T::set_pressure), "", "Pressure")
#ifdef UG_FOR_LUA
			.add_method("set_pressure", static_cast<void (T::*)(const char*)>(&T::set_pressure), "", "Pressure")
#endif

			.add_method("use_elastoplast_mat_behavior", &T::use_elastoplast_mat_behavior)
			.add_method("set_hardening_behavior", &T::set_hardening_behavior)
			.add_method("use_approx_tangent", &T::use_approx_tangent)
			.add_method("init_state_variables", &T::init_state_variables)
			.add_method("stress_eigenvalues_at", &T::stress_eigenvalues_at)
			.add_method("normal_stresses_at", &T::normal_stresses_at)
			.add_method("close_gnuplot_file", &T::close_gnuplot_file)
			.add_method("config_string", &T::config_string)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SmallStrainMechanics", tag);
	}
}

/**
 * Function called for the registration of Dimension dependent parts
 * of the plugin. All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

	//	Material Law Interface
	{
		typedef IMaterialLaw<dim> T;
		string name = string("IMaterialLaw").append(suffix);
		reg.add_class_<T>(name, grp)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "IMaterialLaw", tag);
	}

	//	Hooke Law
	{
		typedef HookeLaw<dim> T;
		typedef IMaterialLaw<dim> TBase;
		string name = string("HookeLaw").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_elasticity_tensor_orthotropic", &T::set_elasticity_tensor_orthotropic, "", "C11#C12#C13#C22#C23#C33#C44#C55#C66")
			.add_method("set_hooke_elasticity_tensor", &T::set_hooke_elasticity_tensor, "", "lambda#mu")
			.add_method("set_hooke_elasticity_tensor_E_nu", &T::set_hooke_elasticity_tensor_E_nu, "", "E#nu")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "HookeLaw", tag);
	}

}


}; // end Functionality

// end group small_strain_mechanics
/// \}

} // end namespace SmallStrainMechanics

/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_SmallStrainMechanics(Registry* reg, string grp)
{
	grp.append("/SpatialDisc/SmallStrainMechanics");
	typedef SmallStrainMechanics::Functionality Functionality;

	try{
		RegisterDimension2d3dDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
