################################################################################\
# SMALL_STRAIN_MECHANICS
#
# The small strain mechanics plugin.
################################################################################\

cmake_minimum_required(VERSION 2.6)

project(SMALL_STRAIN_MECHANICS)

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")


# sources
set(SOURCES	small_strain_mechanics.cpp
			small_strain_mechanics_plugin.cpp)




if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/core/small_strain_mechanics ${SOURCES})
else(buildEmbeddedPlugins)
	add_library(SmallStrainMechanics SHARED ${SOURCES})
	target_link_libraries (SmallStrainMechanics ug4)
endif(buildEmbeddedPlugins)
