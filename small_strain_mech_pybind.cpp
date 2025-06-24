

#include "small_strain_mech_plugin.h"


#ifdef UG_USE_PYBIND11

PYBIND11_MODULE(pysmallstrainmechanics, m)
{
	m.doc() = "SmallStrainMechanics module";
	m.attr("__name__") = "ug4py.smallstrainmechanics";

	ug::pybind::Registry registry(m);
	std::string name("SmallStrainMechanics");

	ug::SmallStrainMechanics::InitUGPlugin(&registry, name);
}
#endif
