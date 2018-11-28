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

#include "small_strain_mech.h"
//#include "adaptive_util.h"
#include "contact/contact.h"

#include "material_laws/hooke.h"
#include "material_laws/damage_law.h"
#include "material_laws/prandtl_reuss.h"
#include "material_laws/mat_law_interface.h"
#include "output_writer/mech_output_writer.h"
#include "bridge/util_overloaded.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace SmallStrainMechanics{


template <typename TDomain>
class DamageFunctionUpdater
{
	public:
		static const int dim = TDomain::dim;
		typedef typename TDomain::grid_type TGrid;
		typedef typename grid_dim_traits<dim>::element_type TElem; 
		typedef typename grid_dim_traits<dim>::side_type TSide; 

		typedef typename TDomain::position_accessor_type TPositionAccessor;

	protected:
		void AveragePositions(	MathVector<dim>& vCenter, 
								const std::vector<MathVector<dim> >& vCornerCoords)
		{
			vCenter = vCornerCoords[0];
			for(size_t j = 1; j < vCornerCoords.size(); ++j)
			{
				vCenter += vCornerCoords[j];
			}
			vCenter *= 1./(number)( vCornerCoords.size());
		}


		void CollectStencilNeighbors(std::vector<TElem*>& vElem,
									 std::vector<DoFIndex>& vIndex,
									 std::vector< MathVector<dim> >& vDistance,
									 TElem* elem,
									 TGrid& grid,
									 TPositionAccessor& aaPos,
									 SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
									 SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0)
		{
			//static const int numNeighborsToFind = 2*dim + (dim * (dim-1)) / 2;
			const size_t fct = 0;

			vElem.clear();
			vIndex.clear();
			vDistance.clear();

			// get vertices of element
			Vertex* const* vVertex = elem->vertices();
			const size_t numVertex = elem->num_vertices();

			// corner coordinates
			std::vector<MathVector<dim> > vCornerCoords;
			vCornerCoords.resize(numVertex);
			for(size_t vrt = 0; vrt < numVertex; ++vrt){
				vCornerCoords[vrt] = aaPos[ vVertex[vrt] ];
			}

			// element midpoint
			MathVector<dim> ElemMidPoint;
			AveragePositions(ElemMidPoint, vCornerCoords);

//			UG_LOG("##############  Element with midpoint " << ElemMidPoint << "    ################   \n")

			// get all sides in order
			typename TGrid::template traits<TSide>::secure_container vSide;
			grid.associated_elements_sorted(vSide, elem);

			//	begin marking
			grid.begin_marking();

			//	mark the elem
			grid.mark(elem);

			////////////////////////////////////////////////////////////////////////////
			// loop all sides
			////////////////////////////////////////////////////////////////////////////
			for(size_t s = 0; s < vSide.size(); ++s){

//				UG_LOG("########  Boundary Side nr " << s << ":   ")

				// get all connected elements
				typename TGrid::template traits<TElem>::secure_container vElemOfSide;
				grid.associated_elements(vElemOfSide, vSide[s]);

				// if no elem found: internal error
				if(vElemOfSide.size() == 0)
					UG_THROW("Huh, why no element?");

				////////////////////////////////////////////////////////////////////////////
				// introduce mirror element if required
				////////////////////////////////////////////////////////////////////////////
				// if only element itself, it must be a boundary element
				if(vElemOfSide.size() == 1){
					vElem.push_back(vElemOfSide[0]);

					// add index
					std::vector<DoFIndex> ind;
					if(spF->inner_dof_indices(elem, fct, ind) != 1) UG_THROW("Wrong number dofs");
					vIndex.push_back(ind[0]);


					// add distance
					vDistance.resize(vDistance.size()+1);
					CollectCornerCoordinates(vCornerCoords, *vSide[s], aaPos);

					MathVector<dim> n;
					ElementNormal<dim>(vSide[s]->reference_object_id(), n, &vCornerCoords[0]);

					ProjectPointToPlane(vDistance.back(), ElemMidPoint, vCornerCoords[0], n);
					VecScaleAdd(vDistance.back(), 2.0, ElemMidPoint, -2.0, vDistance.back());

//					UG_LOG("is boundary, with vDistance: "  << vDistance.back() << "\n");
				}
				////////////////////////////////////////////////////////////////////////////
				// find all direct face neighbors
				////////////////////////////////////////////////////////////////////////////
				// else: add other neighbor and mark
				else{
					for(size_t eos = 0; eos < vElemOfSide.size(); ++eos){

						// neighbor elem
						TElem* neighborElem = vElemOfSide[eos];

						// skip self
						if(neighborElem == elem) continue;

						// add neighbor
						grid.mark(neighborElem);					
						vElem.push_back(neighborElem);

						// add index
						std::vector<DoFIndex> ind;
						if(spF->inner_dof_indices(neighborElem, fct, ind) != 1) UG_THROW("Wrong number dofs");
						vIndex.push_back(ind[0]);

						// add distance
						vDistance.resize(vDistance.size()+1);
						CollectCornerCoordinates(vCornerCoords, *neighborElem, aaPos);
						AveragePositions(vDistance.back(), vCornerCoords);
						VecScaleAdd(vDistance.back(), 1.0, ElemMidPoint, -1.0, vDistance.back());

//						UG_LOG("has neighbor, c"  << vDistance.back() << "\n");

					}
				}
			}

			////////////////////////////////////////////////////////////////////////////
			// add additional elements 
			////////////////////////////////////////////////////////////////////////////
			std::vector<TElem*> vOtherNeighbors;


			int closest = -1;
			number closestDist = std::numeric_limits<number>::max();
			MathVector<dim> distance;

//			UG_LOG("Search for vertices-elems: " << numVertex << "\n");
			for(size_t vrt = 0; vrt < numVertex; ++vrt)
			{
				typename TGrid::template traits<TElem>::secure_container vVertexNeighbor;
				grid.associated_elements(vVertexNeighbor, vVertex[vrt]);

//				UG_LOG(" ++ At vertex "<< vrt << " we have #elems: " << vVertexNeighbor.size() << "\n");
				for(size_t eov = 0; eov < vVertexNeighbor.size(); ++eov)
				{
					TElem* neighborElem = vVertexNeighbor[eov];
					
//					UG_LOG(" ++++ elem "<< eov << " is already marked?: " << grid.is_marked(neighborElem) << "\n");
					if(grid.is_marked(neighborElem)) continue;
					
					grid.mark(neighborElem);
					vOtherNeighbors.push_back(neighborElem);


					CollectCornerCoordinates(vCornerCoords, *neighborElem, aaPos);
					AveragePositions(distance, vCornerCoords);
					VecScaleAppend(distance, -1.0, ElemMidPoint);

					number dist = VecTwoNorm(distance);
//					UG_LOG(" ++++ ++ dist: " << dist << ", distance: "<< distance << ", ElemMidPoint: "<< ElemMidPoint << "\n")
					if(dist < closestDist){
						closest = vOtherNeighbors.size()-1;
					}

				}
			}
			if(dim == 3) UG_THROW("DamageFunctionUpdater: This is 2d only --- extend to 3d by searching for 3 additional neighbors");
			if(closest < 0) UG_THROW("DamageFunctionUpdater: closest not detected.")

			//UG_LOG("vOtherNeighbors.size(): " << vOtherNeighbors.size() << "\n");
			//UG_LOG("closest: " << closest << "\n");
			TElem* otherNeighbor = vOtherNeighbors[closest];
			//UG_LOG("otherNeighbor: " << otherNeighbor << "\n");
			vElem.push_back(otherNeighbor);

			// add index
			std::vector<DoFIndex> ind;
			if(spF->inner_dof_indices(otherNeighbor, fct, ind) != 1) UG_THROW("Wrong number dofs");
			vIndex.push_back(ind[0]);

			// add distance
			vDistance.resize(vDistance.size()+1);
			CollectCornerCoordinates(vCornerCoords, *otherNeighbor, aaPos);
			AveragePositions(vDistance.back(), vCornerCoords);
			VecScaleAdd(vDistance.back(), 1.0, ElemMidPoint, -1.0, vDistance.back());

//			UG_LOG("########  Extra Neighbot with vDistance:  " << vDistance.back() << "\n")

			//	end marking
			grid.end_marking();
		}


	public:
		void init(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0)
		{
			const size_t fct = 0;

			// get domain
			SmartPtr<TDomain> domain = spF->domain();
			SmartPtr<typename TDomain::grid_type> grid = domain->grid();
			typename TDomain::position_accessor_type& aaPos = domain->position_accessor();

			// get dof distribution
			SmartPtr<DoFDistribution> dd = spF->dd();

			//	get iterators
			typename DoFDistribution::traits<TElem>::iterator iter, iterEnd;
			iter = dd->begin<TElem>(SurfaceView::ALL); // SurfaceView::MG_ALL
			iterEnd = dd->end<TElem>(SurfaceView::ALL);

			// resize result vectors
			const size_t numIndex = spF->num_dofs();
			m_vB.resize(numIndex);
			m_vDLambda.resize(numIndex);
			m_vIndex.resize(numIndex);


			// storage (for reuse)
			std::vector<TElem*> vElem;
			std::vector<DoFIndex> vIndex;
			std::vector< MathVector<dim> > vDistance;

			//	loop all vertices
			for(;iter != iterEnd; ++iter)
			{
				//	get vertex
				TElem* elem = *iter;

				////////////////////////////////////////////////////////////////////////////
				// Collect suitable neighbors for stencil creation
				////////////////////////////////////////////////////////////////////////////

				CollectStencilNeighbors(vElem, vIndex, vDistance, elem, *grid, aaPos, spF, spPsi0);

				const int numNeighbors = 2*dim + (dim * (dim-1)) / 2;
				if(vIndex.size() != numNeighbors || vDistance.size() != numNeighbors)
					UG_THROW("Wrong number of neighbors detected: " << vIndex.size());

				////////////////////////////////////////////////////////////////////////////
				// Create interpolation matrix and invert
				////////////////////////////////////////////////////////////////////////////

				DenseMatrix<VariableArray2<number> > BlockInv;
				BlockInv.resize(numNeighbors, numNeighbors);
				BlockInv = 0.0;
				for (size_t j = 0; j < numNeighbors; ++j)
				{
					for (int d = 0; d < dim; ++d)
						BlockInv(j,d) = vDistance[j][d];

					int cnt = dim;
					for (int d1 = 0; d1 < dim-1; ++d1)
						for (int d2 = d1+1; d2 < dim; ++d2)
						{
							BlockInv(j,cnt++) = vDistance[j][d1] * vDistance[j][d2];
						}

					for (int d = 0; d < dim; ++d)
						BlockInv(j,cnt++) = 0.5 * vDistance[j][d] * vDistance[j][d];
				}

				if(!Invert(BlockInv))
					UG_THROW("Cannot invert block");

				////////////////////////////////////////////////////////////////////////////
				// extract second-order derivative subblock
				////////////////////////////////////////////////////////////////////////////


				// add index
				std::vector<DoFIndex> ind;
				if(spF->inner_dof_indices(elem, fct, ind) != 1) UG_THROW("Wrong number dofs");
				const size_t i = ind[0][0];


				m_vB[i].resize(dim, numNeighbors);

				m_vDLambda[i] = 0.0;
				for (int d = 0; d < dim; ++d){
					for (size_t j = 0; j < numNeighbors; ++j){
						m_vB[i](d,j) = BlockInv( (numNeighbors-dim)  + d,  j);

						m_vDLambda[i] -= m_vB[i](d,j);
					}
				}


				m_vIndex[i].resize(vIndex.size());
				for(size_t k = 0; k < vIndex.size(); ++k)
					m_vIndex[i][k] = vIndex[k][0];
			}
		}


		number DLambda(size_t i) {return m_vDLambda[i];}
		number Lambda(size_t i, SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF) {
			const number fi = (*spF)[i];

			number res = 0.0;
			for (size_t j = 0; j < m_vIndex[i].size(); ++j){
				const number fij = (*spF)[ m_vIndex[i][j] ] - fi;

//				UG_LOG("f(" << i << "," << j << "): " << fij << "\n");
				for (int d = 0; d < dim; ++d){
					res += 	fij * m_vB[i](d,j);
//					UG_LOG("m_vB[" << i << "](" << d << "," << j << "): " << m_vB[i](d,j) << "\n");
				}
			}

			return res;
		}

		bool solve(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0,
					const number beta, const number r, const number eps, const int maxIter)
		{
			////////////////////////////////////////////////////////////////////////////
			// check if has to be rebuild 
			////////////////////////////////////////////////////////////////////////////

			// get approximation space
			ConstSmartPtr<ApproximationSpace<TDomain> > approxSpace = spF->approx_space();
			if(approxSpace != spPsi0->approx_space())
				UG_THROW("DamageFunctionUpdater::solve: expected same ApproximationSpace for f and psi0");

			// check revision counter if grid / approx space has changed since last call
			if(m_ApproxSpaceRevision != approxSpace->revision())
			{
				// (re-)initialize setting
				init(spF, spPsi0);

				//	remember revision counter of approx space
				m_ApproxSpaceRevision = approxSpace->revision();
			}



			////////////////////////////////////////////////////////////////////////////
			// apply newton method 
			////////////////////////////////////////////////////////////////////////////

			const size_t numElem = m_vIndex.size();

			number normPhi =  std::numeric_limits<number>::max();
			int iterCnt = 0;
			while(normPhi > eps * numElem && (iterCnt++ <= maxIter) )
			{
				normPhi = 0.0;

				for(size_t i = 0; i < m_vIndex.size(); ++i)
				{
					const number lambda = Lambda(i, spF);
					number& f = (*spF)[i];
					const number& psi0 = (*spPsi0)[i];

					number phi = f * ( psi0  - beta * lambda) - r;

					//if( fabs (lambda - 4.0) > 1e-10  )
					//	UG_LOG("lambda: " << lambda << "\n");
					(*spPsi0)[i] = lambda;
					continue;
//					UG_LOG("DLambda: " << DLambda(i) << "\n");


					if(phi < eps){
//						phi = 0;
//						normPhi += 0.0 * 0.0;
					} else {

/*
						UG_LOG(" ###############  \n");
						UG_LOG("f     : " << f << "\n");
						UG_LOG("psi0  : " << psi0 << "\n");
						UG_LOG("lambda: " << lambda << "\n");
						UG_LOG("DLambda: " << DLambda(i) << "\n");
						UG_LOG("beta  : " << beta << "\n");
						UG_LOG("r     : " << r << "\n");
						UG_LOG("phi   : " << phi << "\n");

*/
						f = f - (phi / (psi0 - beta * (lambda + f * DLambda(i)) ));
//						UG_LOG ("DamageFunctionUpdater: SETTING new value for f: "  << f << "\n");

	
						normPhi += phi*phi;
					}
				}

				normPhi = sqrt(normPhi);
//				UG_LOG(" ######################### (end sweep) ###################################  \n");
//				UG_LOG ("DamageFunctionUpdater: normPhi: "  << normPhi << "\n");	
			}

			if(iterCnt >= maxIter){
				UG_THROW("DamageFunctionUpdater: no convergence after " << iterCnt << " iterations");
			}


			UG_LOG ("DamageFunctionUpdater: normPhi: "  << normPhi << " after " <<iterCnt << " iterations\n");	

			return true;
		}

	protected:

		//	approximation space revision of cached values
		RevisionCounter m_ApproxSpaceRevision;

		std::vector< DenseMatrix<VariableArray2<number> > > m_vB;
		std::vector< number > m_vDLambda;
		std::vector< std::vector<size_t> > m_vIndex;



};




template<typename TDomain>
void MarkDamage(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0,
					IRefiner& refiner,
					number refineFrac, number coarseFrac, 
					number avgRefineFactor, number avgCoarsenFactor,
					int maxLevel)
{
	static const int dim = TDomain::dim;
	typedef typename grid_dim_traits<dim>::element_type TElem; 
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;
	const int fct = 0;

	///////////////////////////
	// statistic
	///////////////////////////

	number maxValue = std::numeric_limits<number>::min();
	number minValue = std::numeric_limits<number>::max();
	number sumValue = 0.0;

	// loop all elems 
	const size_t numElem = spF->size();
	for(size_t i = 0; i < numElem; ++i){
		const number val = (*spF)[i] * (*spPsi0)[i];

		maxValue = std::max(maxValue, val);
		minValue = std::min(minValue, val);
		sumValue += val;
	}
	number avgValue = sumValue / numElem;


	UG_LOG("  +++ numElem: " << numElem << "\n");
	UG_LOG("  +++ maxValue: " << maxValue << "\n");
	UG_LOG("  +++ minValue: " << minValue << "\n");
	UG_LOG("  +++ avgValue: " << avgValue << "\n");


	///////////////////////////
	// selection criteria
	///////////////////////////

	number minValueToRefine = std::min( refineFrac * maxValue, avgRefineFactor * avgValue);
	number maxValueToCoarsen = std::max( (1 + coarseFrac) * minValue, avgCoarsenFactor * avgValue);

	UG_LOG("  +++ minValueToRefine: " << minValueToRefine << "\n");
	UG_LOG("  +++ maxValueToCoarsen: " << maxValueToCoarsen << "\n");

	///////////////////////////
	// marking
	///////////////////////////

	//	reset counter
	int numMarkedRefine = 0, numMarkedCoarse = 0;

	const_iterator iter = spF->template begin<TElem>();
	const_iterator iterEnd = spF->template end<TElem>();

	//	loop elements for marking
	for(; iter != iterEnd; ++iter)
	{
		//	get element
		TElem* elem = *iter;

		std::vector<DoFIndex> ind;
		if(spF->inner_dof_indices(elem, fct, ind) != 1) UG_THROW("Wrong number dofs");
		const number val = DoFRef(*spF, ind[0]) * DoFRef(*spPsi0, ind[0]);

		//	marks for refinement
		if( val >= minValueToRefine)
			if(spF->dd()->multi_grid()->get_level(elem) < maxLevel)
			{
				refiner.mark(elem, RM_REFINE);
				numMarkedRefine++;
			}

		//	marks for coarsening
		if( val <= maxValueToCoarsen)
		{
			refiner.mark(elem, RM_COARSEN);
			numMarkedCoarse++;
		}

	}

	///////////////////////////
	// print infos
	///////////////////////////

#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1)
	{
		UG_LOG("  +++ Marked for refinement on Proc "<<pcl::ProcRank()<<": " << numMarkedRefine << " Elements.\n");
		UG_LOG("  +++ Marked for coarsening on Proc "<<pcl::ProcRank()<<": " << numMarkedCoarse << " Elements.\n");
		pcl::ProcessCommunicator com;
		int numMarkedRefineLocal = numMarkedRefine, numMarkedCoarseLocal = numMarkedCoarse;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
		numMarkedCoarse = com.allreduce(numMarkedCoarseLocal, PCL_RO_SUM);
	}
#endif

	UG_LOG("  +++ Marked for refinement: " << numMarkedRefine << " Elements.\n");
	UG_LOG("  +++ Marked for coarsening: " << numMarkedCoarse << " Elements.\n");
}


template<typename TDomain>
void HadamardProd(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spFPsi0,
					ConstSmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					ConstSmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0)
{
	if(spF->size() != spPsi0->size() || spF->size() != spFPsi0->size())
		UG_THROW("HadamardProd: Size mismatch");

	for(size_t i = 0; i < spF->size(); ++i){
		(*spFPsi0)[i] = (*spF)[i] * (*spPsi0)[i];
	}
}



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
			//.add_method("set_elasticity_tensor_orthotropic2d", &T::set_elasticity_tensor_orthotropic2d, "", "C11#C12#C22#C33")
			.add_method("set_hooke_elasticity_tensor", &T::set_hooke_elasticity_tensor,
					"", "lambda#mu")
			.add_method("set_hooke_elasticity_tensor_E_nu", &T::set_hooke_elasticity_tensor_E_nu,
					"", "E#nu")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "HookeLaw", tag);
	}


	//	Damage Law for Linear Elasticity
	{
		typedef DamageLaw<TDomain> T;
		typedef IMaterialLaw<TDomain> TBase;
		string name = string("DamageLaw").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
//			.add_constructor()
			.template add_constructor<void (*)(SmartPtr<GridFunction<TDomain,CPUAlgebra> >, SmartPtr<GridFunction<TDomain,CPUAlgebra> >)>()
			.add_method("set_elasticity_tensor_orthotropic", &T::set_elasticity_tensor_orthotropic,
					"", "C11#C12#C13#C22#C23#C33#C44#C55#C66")
			//.add_method("set_elasticity_tensor_orthotropic2d", &T::set_elasticity_tensor_orthotropic2d, "", "C11#C12#C22#C33")
			.add_method("set_hooke_elasticity_tensor", &T::set_hooke_elasticity_tensor,
					"", "lambda#mu")
			.add_method("set_hooke_elasticity_tensor_E_nu", &T::set_hooke_elasticity_tensor_E_nu,
					"", "E#nu")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DamageLaw", tag);
	}

	{
		typedef DamageFunctionUpdater<TDomain> T;
		string name = string("DamageFunctionUpdater").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("init", &T::init, "", "")
			.add_method("solve", &T::solve, "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DamageFunctionUpdater", tag);

		reg.add_function("MarkDamage", &MarkDamage<TDomain>, grp);
		reg.add_function("HadamardProd", &HadamardProd<TDomain>, grp);

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
