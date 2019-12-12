/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Raphael Prohl, Andreas Vogel
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

#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/common/geometry_util.h"
#include "lib_grid/refinement/refiner_interface.h"
#include "lib_disc/function_spaces/grid_function_util.h"

#include "small_strain_mech.h"
//#include "adaptive_util.h"
#include "contact/contact.h"

#include "material_laws/hooke.h"
#include "material_laws/scaled_hooke_law.h"
#include "material_laws/skin_law.h"
#include "material_laws/prandtl_reuss.h"
#include "material_laws/mat_law_interface.h"
#include "output_writer/mech_output_writer.h"
#include "bridge/util_overloaded.h"

#include "damage.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace SmallStrainMechanics{

/** 
 *  \defgroup small_strain_mechanics Small Strain Mechanics
 *  \ingroup plugins
 *  This is a plugin for Small Strain Mechanics functionality.
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

	typedef GridFunction<TDomain, TAlgebra> function_type;

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

//	functionality for output
	reg.add_function("normal_stresses_strains",
					&normal_stresses_strains<function_type>, grp);
	reg.add_function("plast_ip",
					&plast_ip<function_type>, grp);
	reg.add_function("equiv_plast_strain",
					&equiv_plast_strain<function_type>, grp);
	//reg.add_function("invariants_kirchhoff_stress",
	//				&invariants_kirchhoff_stress<function_type>, grp);

//	MarkForAdaption_PlasticElem-Indicator
	//reg.add_function("MarkForAdaption_PlasticElem",
	//		 &MarkForAdaption_PlasticElem<TDomain, TAlgebra>, grp);

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
			.add_method("set_material_law", &T::set_material_law, "", "material law")
			.add_method("get_material_law", &T::get_material_law, "", "material law")
			.add_method("set_output_writer", &T::set_output_writer, "", "set output writer")
			.add_method("set_quad_order", &T::set_quad_order, "", "quad order")

			.add_method("set_mass_scale", &T::set_mass_scale, "", "massScale")

			.add_method("set_volume_forces", OVERLOADED_METHOD_PTR(void, T, set_volume_forces, (SmartPtr<CplUserData<MathVector<dim>, dim> >) ) ,"", "Force field")
			.add_method("set_volume_forces", OVERLOADED_METHOD_PTR(void, T, set_volume_forces, (number)), "", "F")
			.add_method("set_volume_forces", OVERLOADED_METHOD_PTR(void, T, set_volume_forces, (number, number)), "", "F_x, F_y")
			.add_method("set_volume_forces", OVERLOADED_METHOD_PTR(void, T, set_volume_forces, (number, number, number)), "", "F_x, F_y, F_z")

			.add_method("set_div_factor", OVERLOADED_METHOD_PTR(void, T, set_div_factor, (SmartPtr<CplUserData<number, dim> >)), "", "Pressure")
			.add_method("set_div_factor", OVERLOADED_METHOD_PTR(void, T, set_div_factor, (number)), "", "Pressure")

			.add_method("set_compress_factor", OVERLOADED_METHOD_PTR(void, T, set_compress_factor, (number)), "", "Compress-Factor")

			.add_method("set_viscous_forces", OVERLOADED_METHOD_PTR(void, T, set_viscous_forces, (SmartPtr<CplUserData<MathVector<dim>, dim> >, SmartPtr<CplUserData<MathVector<dim>, dim> >) ) ,"", "Force field")
#ifdef UG_FOR_LUA
			.add_method("set_volume_forces", OVERLOADED_METHOD_PTR(void, T, set_volume_forces, (const char*)) , "", "Force field")
			.add_method("set_div_factor", OVERLOADED_METHOD_PTR(void, T, set_div_factor, (const char*)), "", "Pressure")
#endif

			.add_method("displacement", &T::displacement)
			.add_method("divergence", &T::divergence)
			.add_method("init_state_variables", &T::init_state_variables)
			.add_method("config_string", &T::config_string)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SmallStrainMechanics", tag);
	}

	//	Material Law Interface
	{
		typedef IMaterialLaw<TDomain> T;
		string name = string("IMaterialLaw").append(suffix);
		reg.add_class_<T>(name, grp)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "IMaterialLaw", tag);
	}

	//	Hooke Law for Linear Elasticity
	{
		typedef HookeLaw<TDomain> T;
		typedef IMaterialLaw<TDomain> TBase;
		string name = string("HookeLaw").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_elasticity_tensor_orthotropic", &T::set_elasticity_tensor_orthotropic,
					"", "C11#C12#C13#C22#C23#C33#C44#C55#C66")
			.add_method("set_elasticity_tensor_orthotropic_E_G_nu", &T::set_elasticity_tensor_orthotropic_E_G_nu, "", "")
			.add_method("set_elasticity_tensor_orthotropic_plain_stress_E_G_nu", &T::set_elasticity_tensor_orthotropic_plain_stress_E_G_nu, "", "")
			.add_method("set_elasticity_tensor_orthotropic_plain_strain_E_G_nu", &T::set_elasticity_tensor_orthotropic_plain_strain_E_G_nu, "", "")
			//.add_method("set_elasticity_tensor_orthotropic2d", &T::set_elasticity_tensor_orthotropic2d, "", "C11#C12#C22#C33")
			.add_method("set_hooke_elasticity_tensor", &T::set_hooke_elasticity_tensor,"", "lambda#mu")
			.add_method("set_hooke_elasticity_tensor_E_nu", &T::set_hooke_elasticity_tensor_E_nu,"", "E#nu")
			.add_method("set_hooke_elasticity_tensor_plain_stress_E_nu", &T::set_hooke_elasticity_tensor_plain_stress_E_nu,"", "E#nu")
			.add_method("set_hooke_elasticity_tensor_plain_strain_E_nu", &T::set_hooke_elasticity_tensor_plain_strain_E_nu,"", "E#nu")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "HookeLaw", tag);
	}

	//	Skin Law for Linear Elasticity
	{
		typedef SkinMaterialLaw<TDomain> T;
		typedef IMaterialLaw<TDomain> TBase;
		string name = string("SkinMaterialLaw").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_hooke_elasticity_tensor", &T::set_hooke_elasticity_tensor,"", "lambda#mu")
			.add_method("set_hooke_elasticity_tensor_E_nu", &T::set_hooke_elasticity_tensor_E_nu,"", "E#nu")
			.add_method("set_hooke_elasticity_tensor_plain_stress_E_nu", &T::set_hooke_elasticity_tensor_plain_stress_E_nu,"", "E#nu")
			.add_method("set_hooke_elasticity_tensor_plain_strain_E_nu", &T::set_hooke_elasticity_tensor_plain_strain_E_nu,"", "E#nu")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SkinMaterialLaw", tag);
	}


	//	Damage Law for Linear Elasticity
	{
		typedef DamageLaw<TDomain> T;
		typedef HookeLaw<TDomain> TBase;
		string name = string("DamageLaw").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<GridFunction<TDomain,CPUAlgebra> >, SmartPtr<GridFunction<TDomain,CPUAlgebra> >)>()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DamageLaw", tag);
	}

	//	TopologyOptimLaw for Linear Elasticity
	{
		typedef TopologyOptimLaw<TDomain> T;
		typedef HookeLaw<TDomain> TBase;
		string name = string("TopologyOptimLaw").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<GridFunction<TDomain,CPUAlgebra> >, SmartPtr<GridFunction<TDomain,CPUAlgebra> >, int)>()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "TopologyOptimLaw", tag);
	}

	{
		typedef DamageFunctionUpdater<TDomain> T;
		string name = string("DamageFunctionUpdater").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("set_disc_type", &T::set_disc_type, "", "")
			.add_method("solve", &T::solve, "", "")
			.add_method("set_quad_rule", &T::set_quad_rule, "", "")
			.add_method("set_debug", &T::set_debug, "", "")
			.add_method("last_num_iterations", &T::last_num_iterations, "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DamageFunctionUpdater", tag);

		reg.add_function("MarkDamage", &MarkDamage<TDomain>, grp);
		reg.add_function("MarkForAdaption_ValueRangeIndicator", &MarkForAdaption_ValueRangeIndicator<TDomain>, grp);
		reg.add_function("HadamardProd", &HadamardProd<TDomain>, grp);
		reg.add_function("MinMaxElementDiameter", &MinMaxElementDiameter<TDomain>, grp);

		reg.add_function("DamageStatistic", &DamageStatistic<TDomain>, grp);
	}


	{
		typedef RelativeDensityUpdater<TDomain> T;
		string name = string("RelativeDensityUpdater").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("set_disc_type", &T::set_disc_type, "", "")
			.add_method("solve", &T::solve, "", "")
			.add_method("set_debug", &T::set_debug, "", "")
			.add_method("set_quad_rule", &T::set_quad_rule, "", "")
			.add_method("set_enforce_local_required_beta", &T::set_enforce_local_required_beta, "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "RelativeDensityUpdater", tag);
	}

	//	Prandtl Reuss Law for small strain ElastoPlasticity
	{
		typedef PrandtlReuss<TDomain> T;
		typedef IMaterialLaw<TDomain> TBase;
		string name = string("PrandtlReuss").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_bulk_modulus", &T::set_bulk_modulus,
					"", "bulkModulus")
			.add_method("set_shear_modulus", &T::set_shear_modulus,
					"", "shearModulus")
			.add_method("set_initial_flow_stress", &T::set_initial_flow_stress,
					"", "initialFlowStress")
			.add_method("set_residual_flow_stress", &T::set_residual_flow_stress,
					"", "residualFlowStress")
			.add_method("set_hardening_modulus", &T::set_hardening_modulus,
					"", "hardeningModulus")
			.add_method("set_hardening_exponent", &T::set_hardening_exponent,
					"", "hardeningExponent")
			.add_method("set_hardening_behavior", &T::set_hardening_behavior,
					"", "hardeningBehavior")
			.add_method("set_tangent_precision", &T::set_tangent_precision,
					"", "tangentPrecision")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "PrandtlReuss", tag);
	}

//	Solid Mechanics Output Writer
   {
		typedef MechOutputWriter<TDomain> T;
		string name = string("MechOutputWriter").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("stress_eigenvalues_at", &T::stress_eigenvalues_at)
			.add_method("normal_stresses_at", &T::normal_stresses_at)
			.add_method("post_timestep", &T::post_timestep)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "MechOutputWriter", tag);
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
/*template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();
}*/

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
		//RegisterDimension2d3dDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
