# Copyright (c) 2012-2014:  G-CSC, Goethe University Frankfurt
# Author: Raphael Prohl
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.




################################################################################\
# SMALL_STRAIN_MECHANICS
#
# The small strain mechanics plugin.
################################################################################\

set(pluginName	SmallStrainMechanics)
set(SOURCES		small_strain_mech.cpp
				small_strain_mech_plugin.cpp)

message(STATUS "# ${pluginName} - BEGIN ##############################################################################")


################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 2.8.12...3.20)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)


################################################################################
# Classic binding (static or dynamic plugin).
################################################################################
if(NOT USE_PYBIND11)
if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} ug4)
endif(buildEmbeddedPlugins)

endif(NOT USE_PYBIND11)

################################################################################
# Python binding (static plugin, dynamic python interface).
################################################################################
if(USE_PYBIND11)

	# include_directories(${SUPERLU_INCLUDE_DIRS})
	
	###############################################
	# 0) Try to locate ug_4 static library....
	###############################################
	set(UG4_STATIC_LIBRARY ug4_s)
	if (BUILD_WITH_PIP)
		get_libug4_from_pip(PYUG4_STATIC_LIBRARY)
		message(STATUS "Library found: ${PYUG4_STATIC_LIBRARY}")
		set(UG4_STATIC_LIBRARY ${PYUG4_STATIC_LIBRARY})
	endif()
	message(STATUS "Using libug4_s: ${UG4_STATIC_LIBRARY}")
	

	###############################################
	# 1) Create (static) library for plugin.
	###############################################
	add_library(${pluginName} STATIC ${SOURCES})

	###############################################
 	# 2) Create (dynamic) python library module.
 	# First argument in 'python_add_library' must match module name in PYBIND11_MODULE call!
	###############################################
	set(PY_TARGET_NAME pysmallstrainmechanics)
	python_add_library(${PY_TARGET_NAME} MODULE small_strain_mech_pybind.cpp ${SOURCES} WITH_SOABI)
	target_link_libraries (${PY_TARGET_NAME} PRIVATE ${pluginName} ${UG4_STATIC_LIBRARY})
	set_target_properties(${PY_TARGET_NAME} PROPERTIES INSTALL_RPATH "$ORIGIN/..:$ORIGIN/../../../lib")
	install(TARGETS ${PY_TARGET_NAME} LIBRARY DESTINATION ug4py COMPONENT pymodules)
endif(USE_PYBIND11)


message(STATUS "# ${pluginName} - END ##############################################################################")
