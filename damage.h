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

	protected:
		void AveragePositions(	MathVector<dim>& vCenter, 
								const std::vector<MathVector<dim> >& vCornerCoords);

		void CollectStencilNeighbors(std::vector<TElem*>& vElem,
									 std::vector<DoFIndex>& vIndex,
									 std::vector< MathVector<dim> >& vDistance,
									 TElem* elem,
									 TGrid& grid,
									 TPositionAccessor& aaPos,
									 SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
									 SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0);


	public:
		void init_ByTaylorExtension(	
				SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
				SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0);

		number DLambda_ByTaylorExtension(size_t i) {return m_vDLambda[i];}
		number Lambda_ByTaylorExtension(size_t i, SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF);


	public:
		void init_ByPartIntegral(	
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0);

		number DLambda_ByPartIntegral(size_t i) {return m_vStencil[i][0];}	
		number Lambda_ByPartIntegral(size_t i, SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF);

	public:
		bool solve(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0,
					const number beta, const number r, 
					const number eps, const int maxIter, const number dampNewton);

		int last_num_iterations() const {return m_lastNumIters;}

		void set_debug(SmartPtr<GridFunctionDebugWriter<TDomain, CPUAlgebra> > spDebugWriter);
		void write_debug(SmartPtr<GridFunction<TDomain, CPUAlgebra> > spGF, std::string name, int call, int iter);
		void write_stencil_matrix_debug(SmartPtr<GridFunction<TDomain, CPUAlgebra> > spGF, std::string name, int call);

	protected:

		SmartPtr<GridFunctionDebugWriter<TDomain, CPUAlgebra> > m_spDebugWriter;

		//	approximation space revision of cached values
		RevisionCounter m_ApproxSpaceRevision;

		std::vector< DenseMatrix<VariableArray2<number> > > m_vB;
		std::vector< number > m_vDLambda;
		std::vector< std::vector<size_t> > m_vIndex;

		std::vector< std::vector<  number > > m_vStencil;

		int m_lastNumIters;
};


////////////////////////////////////////////////////////////////////////////////
// Damage marking
////////////////////////////////////////////////////////////////////////////////


template<typename TDomain>
void MarkDamage(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0,
					IRefiner& refiner,
					number refineFrac, number coarseFrac, 
					number avgRefineFactor, number avgCoarsenFactor,
					int maxLevel);


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
