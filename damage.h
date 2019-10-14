/*
 * Copyright (c) 2019:  Ruhr University Bochum
 * Authors: Andreas Vogel
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

#ifndef __SMALL_STRAIN_MECHANICS__DAMAGE_H_
#define __SMALL_STRAIN_MECHANICS__DAMAGE_H_

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/common/geometry_util.h"
#include "lib_grid/refinement/refiner_interface.h"
#include "lib_disc/function_spaces/grid_function_util.h"

namespace ug{
namespace SmallStrainMechanics{



////////////////////////////////////////////////////////////////////////////////
// Helper traits for constraints
////////////////////////////////////////////////////////////////////////////////

template <int dim> struct contrained_dim_traits;
template <> struct contrained_dim_traits<2>
{
	typedef ConstrainedEdge contrained_side_type;
	typedef ConstrainingEdge contraining_side_type;
};
template <> struct contrained_dim_traits<3>
{
	typedef ConstrainedFace contrained_side_type;
	typedef ConstrainingFace contraining_side_type;
};


template <int dim>
void AveragePositions(	MathVector<dim>& vCenter, 
						const std::vector<MathVector<dim> >& vCornerCoords);


template <typename TDomain>
void CollectSurfaceNeighbors(	
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF, // some dummy function
					typename grid_dim_traits<TDomain::dim>::element_type* elem,
					std::vector< typename grid_dim_traits<TDomain::dim>::element_type* >& vNeighbors);

template <typename TDomain>
void CollectStencilNeighbors_NeumannZeroBND_IndexAndDistance
(
	std::vector< typename grid_dim_traits<TDomain::dim>::element_type* >& vElem,
	std::vector<size_t>& vIndex,
	std::vector< MathVector<TDomain::dim> >& vDistance,
	typename grid_dim_traits<TDomain::dim>::element_type* elem,
	typename TDomain::grid_type& grid,
	typename TDomain::position_accessor_type& aaPos,
	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF
);

template <typename TDomain>
void InitLaplacian_TaylorExpansion(	
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					std::vector< std::vector<  number > >& vStencil,
					std::vector< std::vector<size_t> >& vIndex);

template <typename TDomain>
void InitLaplacian_PartialIntegration(	
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					std::vector< std::vector<  number > >& vStencil,
					std::vector< std::vector<size_t> >& vIndex,
					int quadRuleType, bool fillElemSizeIntoVector = false);


////////////////////////////////////////////////////////////////////////////////
// Damage updater
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
class DamageFunctionUpdater
{
	public:
		static const int dim = TDomain::dim;
		typedef typename TDomain::grid_type TGrid;
		typedef typename grid_dim_traits<dim>::element_type TElem; 
		typedef typename grid_dim_traits<dim>::side_type TSide; 
		typedef typename contrained_dim_traits<dim>::contrained_side_type TContrainedSide; 
		typedef typename contrained_dim_traits<dim>::contraining_side_type TContrainingSide; 
		typedef typename TDomain::position_accessor_type TPositionAccessor;

	/////////////////////////////////////////////////
	// Setup and Solve
	/////////////////////////////////////////////////
	public:
		DamageFunctionUpdater() : m_discType(_PARTIAL_INTEGRATION_), m_quadRuleType(2) {}

		void set_disc_type(const std::string& type);

		void set_quad_rule(int quadRuleType) {m_quadRuleType = quadRuleType;}

		bool solve(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0,
					const number beta, const number r, 
					const number eps, const int maxIter, const number dampNewton);

		int last_num_iterations() const {return m_lastNumIters;}

	protected:
		enum DiscType {_LEAST_SQUARES_, _TAYLOR_EXPANSION_, _PARTIAL_INTEGRATION_};
		int m_discType;
		int m_quadRuleType; // 1 = Midpoint, 2 = Simpson
		int m_lastNumIters;

		RevisionCounter m_ApproxSpaceRevision;   // approximation space revision of cached values

	/////////////////////////////////////////////////
	// Implementation
	/////////////////////////////////////////////////
	protected:
		number DLambda(size_t i) {return m_vStencil[i][0];}	
		number Lambda(size_t i, SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF)
		{
			number res = 0.0;
			for (size_t j = 0; j < m_vIndex[i].size(); ++j)
				res += m_vStencil[i][j] * (*spF)[ m_vIndex[i][j] ];
			return res;
		}

		std::vector< std::vector<  number > > m_vStencil;
		std::vector< std::vector<size_t> > m_vIndex;

	public:
		void set_debug(SmartPtr<GridFunctionDebugWriter<TDomain, CPUAlgebra> > spDebugWriter);

	protected:
		SmartPtr<GridFunctionDebugWriter<TDomain, CPUAlgebra> > m_spDebugWriter;
		void write_debug(SmartPtr<GridFunction<TDomain, CPUAlgebra> > spGF, std::string name, int call, int iter);
		void write_stencil_matrix_debug(SmartPtr<GridFunction<TDomain, CPUAlgebra> > spGF, std::string name, int call);
};


////////////////////////////////////////////////////////////////////////////////
// RelativeDensityUpdater
////////////////////////////////////////////////////////////////////////////////


template <typename TDomain>
class RelativeDensityUpdater
{
	public:
		static const int dim = TDomain::dim;
		typedef typename TDomain::grid_type TGrid;
		typedef typename grid_dim_traits<dim>::element_type TElem; 
		typedef typename grid_dim_traits<dim>::side_type TSide; 
		typedef typename contrained_dim_traits<dim>::contrained_side_type TContrainedSide; 
		typedef typename contrained_dim_traits<dim>::contraining_side_type TContrainingSide; 

		typedef typename TDomain::position_accessor_type TPositionAccessor;

		RelativeDensityUpdater() : m_quadRuleType(2) {}

	public:
		std::vector<number> solve(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spChi,
									SmartPtr<GridFunction<TDomain, CPUAlgebra> > spDrivingForce,
									const number betaStar, const number etaChiStar, 
									const number chiMin, const number dt, const int p,
									const number rho_target, const number MassTol);

		void set_quad_rule(int quadRuleType) {m_quadRuleType = quadRuleType;}

		void set_debug(SmartPtr<GridFunctionDebugWriter<TDomain, CPUAlgebra> > spDebugWriter);

	protected:
		number DLambda(size_t i) {return m_vStencil[i][0];}	
		number Lambda(size_t i, SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF)
		{
			number res = 0.0;
			for (size_t j = 0; j < m_vIndex[i].size(); ++j)
				res += m_vStencil[i][j] * (*spF)[ m_vIndex[i][j] ];
			return res;
		}

	protected:
		SmartPtr<GridFunctionDebugWriter<TDomain, CPUAlgebra> > m_spDebugWriter;
		void write_debug(SmartPtr<GridFunction<TDomain, CPUAlgebra> > spGF, std::string name, int call, int iter);
		void write_stencil_matrix_debug(SmartPtr<GridFunction<TDomain, CPUAlgebra> > spGF, std::string name, int call);

	protected:
		//	approximation space revision of cached values
		RevisionCounter m_ApproxSpaceRevision;

		int m_quadRuleType; // 1 = Midpoint, 2 = Simpson
		std::vector< std::vector<  number > > m_vStencil;
		std::vector< std::vector<size_t> > m_vIndex;


		SmartPtr<GridFunction<TDomain, CPUAlgebra> > m_spElemSize;
		SmartPtr<GridFunction<TDomain, CPUAlgebra> > m_spLaplaceChi;
		SmartPtr<GridFunction<TDomain, CPUAlgebra> > m_spChiTrial;
};


////////////////////////////////////////////////////////////////////////////////
// Damage marking
////////////////////////////////////////////////////////////////////////////////


template<typename TDomain>
void MarkDamage(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0,
					IRefiner& refiner,
					number minValueToRefine, number maxValueToCoarsen,
					int maxLevel,
					const std::vector<MathVector<TDomain::dim,number>* >& vCenter, 
					const std::vector<number>& vRadius);


template<typename TDomain>
std::vector<number> DamageStatistic(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
										SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0);


template<typename TDomain>
void HadamardProd(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spFPsi0,
					ConstSmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					ConstSmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0);


} // end namespace SmallStrainMechanics
}// namespace ug

#include "damage_impl.h"

#endif /* __SMALL_STRAIN_MECHANICS__DAMAGE_H_ */
