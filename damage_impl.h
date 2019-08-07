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

#ifndef __SMALL_STRAIN_MECHANICS__DAMAGE_IMPL_H_
#define __SMALL_STRAIN_MECHANICS__DAMAGE_IMPL_H_

#include "damage.h"
#include "common/math/math_vector_matrix/math_vector_functions.h"

namespace ug{
namespace SmallStrainMechanics{


template <int dim>
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


////////////////////////////////////////////////////////////////////////////////
// By partial integration
////////////////////////////////////////////////////////////////////////////////


template <typename TDomain>
void InitLaplacianByPartialIntegration(	
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF, // some dummy function
					std::vector< std::vector<  number > >& vStencil,
					std::vector< std::vector<size_t> >& vIndex,
					int quadRuleType, bool fillElemSizeIntoVector)
{
	PROFILE_BEGIN_GROUP(InitLaplacianByPartialIntegration, "Small Strain Mech");

	static const int dim = TDomain::dim;
	typedef typename TDomain::grid_type TGrid;
	typedef typename grid_dim_traits<dim>::element_type TElem; 
	typedef typename grid_dim_traits<dim>::side_type TSide; 
	typedef typename contrained_dim_traits<dim>::contrained_side_type TContrainedSide; 
	typedef typename contrained_dim_traits<dim>::contraining_side_type TContrainingSide; 

	typedef typename TDomain::position_accessor_type TPositionAccessor;

	const size_t fct = 0;

	number sideWeight;
	number vertexWeight;

	if(quadRuleType == 1) // Midpoint
	{
		sideWeight = 1.0;
		vertexWeight = 0.0;
	} 
	else if (quadRuleType == 2) // Simpson's
	{
		sideWeight = 4.0/6.0;
		vertexWeight = 1.0/6.0;		
	} 
	else
		UG_THROW("InitLaplacianByPartialIntegration: wrong quad rule")


	// get domain
	SmartPtr<TDomain> domain = spF->domain();
	SmartPtr<typename TDomain::grid_type> spGrid = domain->grid();
	typename TDomain::grid_type& grid = *spGrid;
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();

	// get dof distribution
	SmartPtr<DoFDistribution> dd = spF->dd();

	//	get iterators
	typename DoFDistribution::traits<TElem>::iterator iter, iterEnd;
	iter = dd->begin<TElem>(SurfaceView::ALL); // SurfaceView::MG_ALL
	iterEnd = dd->end<TElem>(SurfaceView::ALL);

	// resize result vectors
	const size_t numIndex = spF->num_dofs();
	vIndex.resize(numIndex);
	vStencil.resize(numIndex);

	// storage (for reuse)
	std::vector<MathVector<dim> > vCornerCoords, vNbrCornerCoords;
	MathVector<dim> ElemMidPoint, Normal, Distance;

	//	loop all vertices
	for(;iter != iterEnd; ++iter)
	{
		//	get vertex
		TElem* elem = *iter;

		// store index
		std::vector<DoFIndex> ind;
		if(spF->inner_dof_indices(elem, fct, ind) != 1) UG_THROW("Wrong number dofs");
		const size_t i = ind[0][0];

		// add self-coupling
		vIndex[i].clear();
		vIndex[i].push_back(i);				

		// reset stencil
		vStencil[i].clear(); 
		vStencil[i].resize(1, 0.0);

		// get vertices of element
		Vertex* const* vVertex = elem->vertices();
		const size_t numVertex = elem->num_vertices();

		// corner coordinates
		vCornerCoords.resize(numVertex);
		for(size_t vrt = 0; vrt < numVertex; ++vrt){
			vCornerCoords[vrt] = aaPos[ vVertex[vrt] ];
		}

		// element volume
		const number vol = ElementSize<dim>(elem->reference_object_id(), &vCornerCoords[0]);

		if(fillElemSizeIntoVector)
			(*spF)[i] = vol;

		// element midpoint
		AveragePositions<dim>(ElemMidPoint, vCornerCoords);

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

			// side normal
			SideNormal<dim>(elem->reference_object_id(), Normal, s, &vCornerCoords[0]);

			// get all connected elements
			typename TGrid::template traits<TElem>::secure_container vElemOfSide;
			grid.associated_elements(vElemOfSide, vSide[s]);

			// if no elem found: internal error
			if(vElemOfSide.size() == 0)
				UG_THROW("Huh, why no element? Should be at least the considered elem itself");

			////////////////////////////////////////////////////////////////////////////
			// handle sides with only one element
			// ( if only one element for side, it must be a boundary element or contrained)
			////////////////////////////////////////////////////////////////////////////
			if(vElemOfSide.size() == 1){

				////////////////////////////////////////////////////////////////////////////
				// handle constraint edge
				// 
				////////////////////////////////////////////////////////////////////////////
				if(vSide[s]->is_constrained()){

					TContrainedSide* cSide = dynamic_cast<TContrainedSide*>(vSide[s]);					
					TSide* constrainingSide = dynamic_cast<TSide*>(cSide->get_constraining_object());

					typename TGrid::template traits<TElem>::secure_container vElemOfContrainingSide;
					grid.associated_elements(vElemOfContrainingSide, constrainingSide);
					if(vElemOfContrainingSide.size() != 2) UG_THROW("Huh, should be 2 at constraining side");

					for(size_t nbr = 0; nbr < vElemOfContrainingSide.size(); ++nbr){

						TElem* neighborElem = vElemOfContrainingSide[nbr];
						if(grid.template num_children<TElem,TElem>(neighborElem) > 0) continue;

						grid.mark(neighborElem);					

						// add distance
						CollectCornerCoordinates(vNbrCornerCoords, *neighborElem, aaPos);
						AveragePositions<dim>(Distance, vNbrCornerCoords);
						VecScaleAdd(Distance, 1.0, ElemMidPoint, -1.0, Distance);

						const number cpl = (sideWeight) * VecDot(Normal, Distance) / ( VecTwoNormSq(Distance) * vol);
						vStencil[i].push_back( -cpl );
						vStencil[i][0]   +=  cpl;

						std::vector<DoFIndex> ind;
						if(spF->inner_dof_indices(neighborElem, fct, ind) != 1) UG_THROW("Wrong number dofs");
						vIndex[i].push_back(ind[0][0]);

					}
				} 
				////////////////////////////////////////////////////////////////////////////
				// handle mirror elements
				// ( if only one element for side, it must be a boundary element)
				////////////////////////////////////////////////////////////////////////////
				else {

					ProjectPointToPlane(Distance, ElemMidPoint, aaPos[ (vSide[s]->vertex(0)) ], Normal);
					VecScaleAdd(Distance, 2.0, ElemMidPoint, -2.0, Distance);

					const number cpl = VecDot(Normal, Distance) / ( VecTwoNormSq(Distance) * vol);
					vStencil[i].push_back( -cpl );
					vStencil[i][0]   +=  cpl;

					std::vector<DoFIndex> ind;
					if(spF->inner_dof_indices(elem, fct, ind) != 1) UG_THROW("Wrong number dofs");
					vIndex[i].push_back(ind[0][0]);
				}
			}
			////////////////////////////////////////////////////////////////////////////
			// handle hanging sides with 2 elements
			// ( coupled to all fine neighbor elems )
			////////////////////////////////////////////////////////////////////////////
			else if(vSide[s]->is_constraining()){

				TContrainingSide* cSide = dynamic_cast<TContrainingSide*>(vSide[s]);
				const size_t numConstrained = cSide->template num_constrained<TSide>();

				for(size_t cs = 0; cs < numConstrained; ++cs){

					TSide* constrainedSide = cSide->template constrained<TSide>(cs);

					// neighbor elem
					typename TGrid::template traits<TElem>::secure_container vElemOfContrainedSide;
					grid.associated_elements(vElemOfContrainedSide, constrainedSide);
					if(vElemOfContrainedSide.size() != 1) UG_THROW("Huh, should be 1 at constrained side");
					TElem* neighborElem = vElemOfContrainedSide[0];
					grid.mark(neighborElem);					

					// add distance
					CollectCornerCoordinates(vNbrCornerCoords, *neighborElem, aaPos);
					AveragePositions<dim>(Distance, vNbrCornerCoords);
					VecScaleAdd(Distance, 1.0, ElemMidPoint, -1.0, Distance);

					const number cpl = (sideWeight / numConstrained) * VecDot(Normal, Distance) / ( VecTwoNormSq(Distance) * vol);
					vStencil[i].push_back( -cpl );
					vStencil[i][0]   +=  cpl;

					std::vector<DoFIndex> ind;
					if(spF->inner_dof_indices(neighborElem, fct, ind) != 1) UG_THROW("Wrong number dofs");
					vIndex[i].push_back(ind[0][0]);
				}

			}
			////////////////////////////////////////////////////////////////////////////
			// regular refined sides
			// (find all direct face neighbors)
			////////////////////////////////////////////////////////////////////////////
			else{
				for(size_t eos = 0; eos < vElemOfSide.size(); ++eos){

					// if more than 2 elem found: internal error
					if(vElemOfSide.size() != 2)
						UG_THROW("Huh, why more than 2 elements of side?");

					// neighbor elem
					TElem* neighborElem = vElemOfSide[eos];

					// check
					if(grid.template num_children<TElem,TElem>(neighborElem) > 0) 
						UG_THROW("Huh, why not on top level?");

					// skip self
					if(neighborElem == elem) continue;
					grid.mark(neighborElem);					

					// add distance
					CollectCornerCoordinates(vNbrCornerCoords, *neighborElem, aaPos);
					AveragePositions<dim>(Distance, vNbrCornerCoords);
					VecScaleAdd(Distance, 1.0, ElemMidPoint, -1.0, Distance);

					const number cpl = sideWeight * VecDot(Normal, Distance) / ( VecTwoNormSq(Distance) * vol);
					vStencil[i].push_back( -cpl );
					vStencil[i][0]   +=  cpl;

					std::vector<DoFIndex> ind;
					if(spF->inner_dof_indices(neighborElem, fct, ind) != 1) UG_THROW("Wrong number dofs");
					vIndex[i].push_back(ind[0][0]);
				}
			}
		
		} // end loop sides

		////////////////////////////////////////////////////////////////////////////
		// loop all sides (FOR DIAG COUPLING)
		////////////////////////////////////////////////////////////////////////////
		if(vertexWeight != 0.0){

			for(size_t s = 0; s < vSide.size(); ++s){

				// side normal
				SideNormal<dim>(elem->reference_object_id(), Normal, s, &vCornerCoords[0]);

				// get all connected elements
				typename TGrid::template traits<TElem>::secure_container vElemOfSide;
				grid.associated_elements(vElemOfSide, vSide[s]);

				// if no elem found: internal error
				if(vElemOfSide.size() == 0)
					UG_THROW("Huh, why no element? Should be at least the considered elem itself");


				////////////////////////////////////////////////////////////////////////////
				// diagonal elements (only coupled via a common vertex)
				// (find all vertex coupled neighbors)
				////////////////////////////////////////////////////////////////////////////
				
				// loop vertices of side
				Vertex* const* vVertex = vSide[s]->vertices();
				const size_t numVertex = vSide[s]->num_vertices();


				for(size_t vrt = 0; vrt < numVertex; ++vrt)
				{
					std::vector<TElem*> vDiagNeighbors;

					// collect all diag leaf-elems on same level
					typename TGrid::template traits<TElem>::secure_container vVertexNeighbor;
					grid.associated_elements(vVertexNeighbor, vVertex[vrt]);
					for(size_t eov = 0; eov < vVertexNeighbor.size(); ++eov)
					{
						TElem* neighborElem = vVertexNeighbor[eov];
						if(grid.is_marked(neighborElem)) continue;			
						if(grid.template num_children<TElem,TElem>(neighborElem) > 0) continue;

						//UG_LOG("Over Same level\n")
			
						//grid.mark(neighborElem);

						vDiagNeighbors.push_back(neighborElem);
					}

					// collect all diag leaf-elems on finer level
					if(grid.template num_children<Vertex,Vertex>(vVertex[vrt]) > 0){

						if(grid.template num_children<Vertex,Vertex>(vVertex[vrt]) != 1)
							UG_THROW("Huh, why more than one vertex child?")

						Vertex* fineVrt = grid.template get_child<Vertex,Vertex>(vVertex[vrt], 0);

						typename TGrid::template traits<TElem>::secure_container vVertexNeighbor;
						grid.associated_elements(vVertexNeighbor, fineVrt);
						for(size_t eov = 0; eov < vVertexNeighbor.size(); ++eov)
						{
							TElem* neighborElem = vVertexNeighbor[eov];
							if(grid.is_marked(neighborElem)) continue;			
							if(grid.template num_children<TElem,TElem>(neighborElem) > 0) continue;

							//UG_LOG("Over Children\n")
				
							//grid.mark(neighborElem);

							vDiagNeighbors.push_back(neighborElem);
						}

					}

					// collect all diag leaf-elems on coarser level
					if(grid.get_parent(vVertex[vrt]) != 0){


						Vertex* coarseVrt = dynamic_cast<Vertex*>(grid.get_parent(vVertex[vrt]));

						if(coarseVrt != 0){

							typename TGrid::template traits<TElem>::secure_container vVertexNeighbor;
							grid.associated_elements(vVertexNeighbor, coarseVrt);
							for(size_t eov = 0; eov < vVertexNeighbor.size(); ++eov)
							{
								TElem* neighborElem = vVertexNeighbor[eov];
								if(grid.is_marked(neighborElem)) continue;			
								if(grid.template num_children<TElem,TElem>(neighborElem) > 0) continue;
					
								//UG_LOG("Over Parent\n")

								//grid.mark(neighborElem);

								vDiagNeighbors.push_back(neighborElem);
							}
						}
					}

					// add contribution
					const size_t numDiagElem = vDiagNeighbors.size();
					for(size_t diag = 0; diag < numDiagElem; ++diag){

						TElem* neighborElem = vDiagNeighbors[diag];

						// add distance
						CollectCornerCoordinates(vNbrCornerCoords, *neighborElem, aaPos);
						AveragePositions<dim>(Distance, vNbrCornerCoords);
						VecScaleAdd(Distance, 1.0, ElemMidPoint, -1.0, Distance);


						//UG_LOG("Elem: "<<i<<", side: "<<s<<", vrt: "<<vrt<<", Distance: "<<Distance<<", Normal: "<<Normal<<"\n");

						if(dim == 3)
							UG_THROW("This implementation is 2d only, currently. Handle vertex neighbors properly in 3d...");

						std::vector<DoFIndex> ind;
						if(spF->inner_dof_indices(neighborElem, fct, ind) != 1) UG_THROW("Wrong number dofs");
						const size_t j = ind[0][0];

						size_t k = 0;
						for(; k < vIndex[i].size(); ++k){
							if(vIndex[i][k] == j) 
								break;
						}
						if(k == vIndex[i].size()) 
						{
							vIndex[i].push_back(j);
							vStencil[i].push_back(0.0);
							// k == vIndex[i].size();
						}

						const number cpl = (vertexWeight / numDiagElem) * VecDot(Normal, Distance) / ( VecTwoNormSq(Distance) * vol);
						vStencil[i][k]   -= cpl;
						vStencil[i][0]   +=  cpl;
					}
				}
			} // end side loop

		} // end diag-elem block

		//	end marking
		grid.end_marking();

	} // end element loop

	static int call = 0; call++;
	//write_stencil_matrix_debug(spF, "Stencil", call);

}






////////////////////////////////////////////////////////////////////////////////
// Collect Surface Neighbors
////////////////////////////////////////////////////////////////////////////////



template <typename TDomain>
void CollectSurfaceNeighbors(	
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF, // some dummy function
					typename grid_dim_traits<TDomain::dim>::element_type* elem,
					std::vector< typename grid_dim_traits<TDomain::dim>::element_type * >& vNeighbors)
{
	PROFILE_BEGIN_GROUP(InitLaplacianByPartialIntegration, "Small Strain Mech");

	static const int dim = TDomain::dim;
	typedef typename TDomain::grid_type TGrid;
	typedef typename grid_dim_traits<dim>::element_type TElem; 
	typedef typename grid_dim_traits<dim>::side_type TSide; 
	typedef typename contrained_dim_traits<dim>::contrained_side_type TContrainedSide; 
	typedef typename contrained_dim_traits<dim>::contraining_side_type TContrainingSide; 

	typedef typename TDomain::position_accessor_type TPositionAccessor;

	// get domain
	SmartPtr<TDomain> domain = spF->domain();
	SmartPtr<typename TDomain::grid_type> spGrid = domain->grid();
	typename TDomain::grid_type& grid = *spGrid;

	// get dof distribution
	SmartPtr<DoFDistribution> dd = spF->dd();

	//	get iterators
	typename DoFDistribution::traits<TElem>::iterator iter, iterEnd;
	iter = dd->begin<TElem>(SurfaceView::ALL); // SurfaceView::MG_ALL
	iterEnd = dd->end<TElem>(SurfaceView::ALL);


	// clear container
	vNeighbors.clear();

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

		// get all connected elements
		typename TGrid::template traits<TElem>::secure_container vElemOfSide;
		grid.associated_elements(vElemOfSide, vSide[s]);

		// if no elem found: internal error
		if(vElemOfSide.size() == 0)
			UG_THROW("Huh, why no element? Should be at least the considered elem itself");

		////////////////////////////////////////////////////////////////////////////
		// handle sides with only one element
		// ( if only one element for side, it must be a boundary element or contrained)
		////////////////////////////////////////////////////////////////////////////
		if(vElemOfSide.size() == 1){

			////////////////////////////////////////////////////////////////////////////
			// handle constraint edge
			// 
			////////////////////////////////////////////////////////////////////////////
			if(vSide[s]->is_constrained()){

				TContrainedSide* cSide = dynamic_cast<TContrainedSide*>(vSide[s]);					
				TSide* constrainingSide = dynamic_cast<TSide*>(cSide->get_constraining_object());

				typename TGrid::template traits<TElem>::secure_container vElemOfContrainingSide;
				grid.associated_elements(vElemOfContrainingSide, constrainingSide);
				if(vElemOfContrainingSide.size() != 2) UG_THROW("Huh, should be 2 at constraining side");

				for(size_t nbr = 0; nbr < vElemOfContrainingSide.size(); ++nbr){

					TElem* neighborElem = vElemOfContrainingSide[nbr];
					if(grid.template num_children<TElem,TElem>(neighborElem) > 0) continue;

					grid.mark(neighborElem);					
					vNeighbors.push_back(neighborElem);
				} 
			}
		}
		////////////////////////////////////////////////////////////////////////////
		// handle hanging sides with 2 elements
		// ( coupled to all fine neighbor elems )
		////////////////////////////////////////////////////////////////////////////
		else if(vSide[s]->is_constraining()){

			TContrainingSide* cSide = dynamic_cast<TContrainingSide*>(vSide[s]);
			const size_t numConstrained = cSide->template num_constrained<TSide>();

			for(size_t cs = 0; cs < numConstrained; ++cs){

				TSide* constrainedSide = cSide->template constrained<TSide>(cs);

				// neighbor elem
				typename TGrid::template traits<TElem>::secure_container vElemOfContrainedSide;
				grid.associated_elements(vElemOfContrainedSide, constrainedSide);
				if(vElemOfContrainedSide.size() != 1) UG_THROW("Huh, should be 1 at constrained side");
				TElem* neighborElem = vElemOfContrainedSide[0];
				grid.mark(neighborElem);					
				vNeighbors.push_back(neighborElem);

			}

		}
		////////////////////////////////////////////////////////////////////////////
		// regular refined sides
		// (find all direct face neighbors)
		////////////////////////////////////////////////////////////////////////////
		else{
			for(size_t eos = 0; eos < vElemOfSide.size(); ++eos){

				// if more than 2 elem found: internal error
				if(vElemOfSide.size() != 2)
					UG_THROW("Huh, why more than 2 elements of side?");

				// neighbor elem
				TElem* neighborElem = vElemOfSide[eos];

				// check
				if(grid.template num_children<TElem,TElem>(neighborElem) > 0) 
					UG_THROW("Huh, why not on top level?");

				// skip self
				if(neighborElem == elem) continue;
				grid.mark(neighborElem);					
				vNeighbors.push_back(neighborElem);
			}
		}
		
	} // end loop sides

	////////////////////////////////////////////////////////////////////////////
	// loop all sides (FOR DIAG COUPLING)
	////////////////////////////////////////////////////////////////////////////
	for(size_t s = 0; s < vSide.size(); ++s){

		// get all connected elements
		typename TGrid::template traits<TElem>::secure_container vElemOfSide;
		grid.associated_elements(vElemOfSide, vSide[s]);

		// if no elem found: internal error
		if(vElemOfSide.size() == 0)
			UG_THROW("Huh, why no element? Should be at least the considered elem itself");


		////////////////////////////////////////////////////////////////////////////
		// diagonal elements (only coupled via a common vertex)
		// (find all vertex coupled neighbors)
		////////////////////////////////////////////////////////////////////////////
		
		// loop vertices of side
		Vertex* const* vVertex = vSide[s]->vertices();
		const size_t numVertex = vSide[s]->num_vertices();


		for(size_t vrt = 0; vrt < numVertex; ++vrt)
		{
			// collect all diag leaf-elems on same level
			typename TGrid::template traits<TElem>::secure_container vVertexNeighbor;
			grid.associated_elements(vVertexNeighbor, vVertex[vrt]);
			for(size_t eov = 0; eov < vVertexNeighbor.size(); ++eov)
			{
				TElem* neighborElem = vVertexNeighbor[eov];
				if(grid.is_marked(neighborElem)) continue;			
				if(grid.template num_children<TElem,TElem>(neighborElem) > 0) continue;

				grid.mark(neighborElem);
				vNeighbors.push_back(neighborElem);
			}

			// collect all diag leaf-elems on finer level
			if(grid.template num_children<Vertex,Vertex>(vVertex[vrt]) > 0){

				if(grid.template num_children<Vertex,Vertex>(vVertex[vrt]) != 1)
					UG_THROW("Huh, why more than one vertex child?")

				Vertex* fineVrt = grid.template get_child<Vertex,Vertex>(vVertex[vrt], 0);

				typename TGrid::template traits<TElem>::secure_container vVertexNeighbor;
				grid.associated_elements(vVertexNeighbor, fineVrt);
				for(size_t eov = 0; eov < vVertexNeighbor.size(); ++eov)
				{
					TElem* neighborElem = vVertexNeighbor[eov];
					if(grid.is_marked(neighborElem)) continue;			
					if(grid.template num_children<TElem,TElem>(neighborElem) > 0) continue;
		
					grid.mark(neighborElem);
					vNeighbors.push_back(neighborElem);
				}

			}

			// collect all diag leaf-elems on coarser level
			if(grid.get_parent(vVertex[vrt]) != 0){


				Vertex* coarseVrt = dynamic_cast<Vertex*>(grid.get_parent(vVertex[vrt]));

				if(coarseVrt != 0){

					typename TGrid::template traits<TElem>::secure_container vVertexNeighbor;
					grid.associated_elements(vVertexNeighbor, coarseVrt);
					for(size_t eov = 0; eov < vVertexNeighbor.size(); ++eov)
					{
						TElem* neighborElem = vVertexNeighbor[eov];
						if(grid.is_marked(neighborElem)) continue;			
						if(grid.template num_children<TElem,TElem>(neighborElem) > 0) continue;
			
						grid.mark(neighborElem);
						vNeighbors.push_back(neighborElem);
					}
				}
			}


		}
	} // end side loop


	//	end marking
	grid.end_marking();

}

////////////////////////////////////////////////////////////////////////////////
// By taylor expansion
////////////////////////////////////////////////////////////////////////////////
/*

template <typename TDomain>
void DamageFunctionUpdater<TDomain>::
CollectStencilNeighbors(std::vector<TElem*>& vElem,
						 std::vector<DoFIndex>& vIndex,
						 std::vector< MathVector<dim> >& vDistance,
						 TElem* elem,
						 TGrid& grid,
						 TPositionAccessor& aaPos,
						 SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
						 SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0)
{
	PROFILE_BEGIN_GROUP(DamageFunctionUpdater_CollectStencilNeighbors, "Small Strain Mech");

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
	AveragePositions<dim>(ElemMidPoint, vCornerCoords);

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

				// if more than 2 elem found: internal error
				if(vElemOfSide.size() != 2)
					UG_THROW("Huh, why more than 2 elements of side?");

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
				AveragePositions<dim>(vDistance.back(), vCornerCoords);
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
			AveragePositions<dim>(distance, vCornerCoords);
			VecScaleAppend(distance, -1.0, ElemMidPoint);

			number dist = VecTwoNorm(distance);
//					UG_LOG(" ++++ ++ dist: " << dist << ", distance: "<< distance << ", ElemMidPoint: "<< ElemMidPoint << "\n")
			if(dist < closestDist){
				closest = vOtherNeighbors.size()-1;
				closestDist = dist;
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
	AveragePositions<dim>(vDistance.back(), vCornerCoords);
	VecScaleAdd(vDistance.back(), 1.0, ElemMidPoint, -1.0, vDistance.back());

//			UG_LOG("########  Extra Neighbot with vDistance:  " << vDistance.back() << "\n")

	//	end marking
	grid.end_marking();
}


template <typename TDomain>
void DamageFunctionUpdater<TDomain>::
init_ByTaylorExtension(	
				SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
				SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0)
{
	PROFILE_BEGIN_GROUP(DamageFunctionUpdater_init_ByTaylorExtension, "Small Strain Mech");

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


///////////// BEBUG (begin) ///////////////
	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spCond_F = spF->clone();
	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spCond_1 = spF->clone();
	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spCond_infty = spF->clone();
///////////// BEBUG (end) ///////////////

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

			if(cnt != numNeighbors)
				UG_THROW("Wrong number of equations")
		}


		DenseMatrix<VariableArray2<number> > Block;
		Block = BlockInv;



		if(!Invert(BlockInv))
			UG_THROW("Cannot invert block");


///////////// BEBUG (begin) ///////////////
		for (size_t i = 0; i < numNeighbors; ++i){
			for (size_t j = 0; j < numNeighbors; ++j){

				number res = 0.0;
				for(size_t k = 0; k < numNeighbors; ++k)
					res += Block(i,k) * BlockInv(k,j);

				if(i == j){
					if( fabs(res - 1.0) > 1e-10 )
						UG_THROW("Inverse wrong, should be 1.0, but is: " << res);
				} else {
					if( fabs(res - 0.0) > 1e-10 )
						UG_THROW("Inverse wrong, should be 0.0, but is: " << res);							
				}
			}
		}

///////////// BEBUG (end) ///////////////



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

//						if(vIndex[j][0] != i){
					m_vB[i](d,j) = BlockInv( (numNeighbors-dim)  + d,  j);
					m_vDLambda[i] -= m_vB[i](d,j);
//						} else {
//							m_vB[i](d,j) = 0.0; 
//						}

			}
		}


///////////// BEBUG (begin) ///////////////
		const number D_F = MatFrobeniusNorm(Block);
		const number DInv_F = MatFrobeniusNorm(BlockInv);

		const number D_1 = MatOneNorm(Block);
		const number DInv_1 = MatOneNorm(BlockInv);

		const number D_infty = MatInftyNorm(Block);
		const number DInv_infty = MatInftyNorm(BlockInv);

		(*spCond_F)[i] = D_F * DInv_F;
		(*spCond_1)[i] = D_1 * DInv_1;
		(*spCond_infty)[i] = D_infty * DInv_infty;
///////////// BEBUG (end) ///////////////

//				UG_LOG(" ++++++ Distances  Elem "  << i << "\n")
//		for (int k = 0; k < numNeighbors; ++k)
//			UG_LOG(" ++"  << k << ": "  << vDistance[k] <<"\n")
//		UG_LOG(" +++ DLambda "  << m_vDLambda[i] << "\n")
//
		
		m_vIndex[i].resize(vIndex.size());
		for(size_t k = 0; k < vIndex.size(); ++k)
			m_vIndex[i][k] = vIndex[k][0];
	}


///////////// BEBUG (begin) ///////////////
	write_debug(spCond_F, "Cond_D_F", 0, 0);
	write_debug(spCond_1, "Cond_D_1", 0, 0);
	write_debug(spCond_infty, "Cond_D_infty", 0, 0);
///////////// BEBUG (end) ///////////////

}

template <typename TDomain>
number DamageFunctionUpdater<TDomain>::
Lambda_ByTaylorExtension(size_t i, SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF) 
{
	const number fi = (*spF)[i];

	number res = 0.0;
	for (size_t j = 0; j < m_vIndex[i].size(); ++j){
		const number fij = (*spF)[ m_vIndex[i][j] ] - fi;

		for (int d = 0; d < dim; ++d){
			res += 	fij * m_vB[i](d,j);
		}
	}
	return res;
}
*/



////////////////////////////////////////////////////////////////////////////////
// DamageFunctionUpdater
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
bool DamageFunctionUpdater<TDomain>::
solve(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
		SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0,
		const number beta, const number r, 
		const number eps, const int maxIter, const number dampNewton)
{
	PROFILE_BEGIN_GROUP(DamageFunctionUpdater_solve, "Small Strain Mech");

	////////////////////////////////////////////////////////////////////////////
	// check if has to be rebuild 
	////////////////////////////////////////////////////////////////////////////

	static int call = 0; call++;	

	// get approximation space
	ConstSmartPtr<ApproximationSpace<TDomain> > approxSpace = spF->approx_space();
	if(approxSpace != spPsi0->approx_space())
		UG_THROW("DamageFunctionUpdater<TDomain>::solve: expected same ApproximationSpace for f and psi0");

	// check revision counter if grid / approx space has changed since last call
	if(m_ApproxSpaceRevision != approxSpace->revision())
	{
		// (re-)initialize setting
		InitLaplacianByPartialIntegration(spF, m_vStencil, m_vIndex, m_quadRuleType);

		//	remember revision counter of approx space
		m_ApproxSpaceRevision = approxSpace->revision();
	}



	////////////////////////////////////////////////////////////////////////////
	// apply newton method 
	////////////////////////////////////////////////////////////////////////////

//			const size_t numElem = m_vIndex.size();
//			const number sqrtNumElem = sqrt(numElem);

	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spLambdaOld = spF->clone();
	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPhi = spF->clone();


	// number normPhi =  std::numeric_limits<number>::max();
	int iterCnt = 0;

///////////// BEBUG (begin) ///////////////
/*	for(size_t i = 0; i < m_vIndex.size(); ++i)
		(*spLambdaOld)[i] = DLambda(i);
	write_debug(spLambdaOld, "DLambda", call, iterCnt);

	write_debug(spPsi0, "Psi0", call, iterCnt);
	write_debug(spF, "F", call, iterCnt);
*///////////// BEBUG (end) ///////////////

	number maxPhi = std::numeric_limits<number>::max();

//	while(normPhi > eps * sqrtNumElem && (iterCnt++ <= maxIter) )
	while(maxPhi > eps && (iterCnt++ <= maxIter) )
	{
		// normPhi = 0.0;
		maxPhi  = 0.0;

		for(size_t i = 0; i < m_vIndex.size(); ++i)
			(*spLambdaOld)[i] = Lambda(i, spF);

///////////// BEBUG (begin) ///////////////
//		write_debug(spLambdaOld, "Lambda", call, iterCnt);
///////////// BEBUG (end) ///////////////

		for(size_t i = 0; i < m_vIndex.size(); ++i)
		{
			const number lambda = (*spLambdaOld)[i];
			const number& psi0 = (*spPsi0)[i];
			number& f = (*spF)[i];

			number phi = f * ( psi0  - beta * lambda) - r;

///////////// BEBUG (begin) ///////////////
//			(*spPhi)[i] = phi;
///////////// BEBUG (end) ///////////////


			if(phi < eps){
//						phi = 0;
//						normPhi += 0.0 * 0.0;
			} else {

				// ORIGINAL
				// f = f - (phi / (psi0 - beta * (lambda + f * DLambda(i)) ));


				// DAMPED NEWTON
				f = f - dampNewton * (phi / (psi0 - beta * (lambda + f * DLambda(i)) ));

				// QUASI-NEWTON
				//f = f - (1/10)*(phi / (psi0 - beta * lambda));

				//  FIXPOINT
				//f =  phiScale * (f * ( psi0  - beta * lambda) - r) + f;


///////////// BEBUG (begin) ///////////////
				if(std::isfinite(f) == false){

					UG_LOG(" ###############  \n");
					UG_LOG("f     : " << f << "\n");
					UG_LOG("psi0  : " << psi0 << "\n");
					UG_LOG("lambda: " << lambda << "\n");
					UG_LOG("DLambda: " << DLambda(i) << "\n");
					UG_LOG("beta  : " << beta << "\n");
					UG_LOG("r     : " << r << "\n");
					UG_LOG("phi   : " << phi << "\n");
					UG_LOG(" ###############  \n");
					UG_THROW("Value for f not finite, but: " << f);
				}
////////////// BEBUG (end) ///////////////

				// normPhi += phi*phi;
				maxPhi = std::max(maxPhi, phi);
			}
		}

///////////// BEBUG (begin) ///////////////
//		write_debug(spPhi, "Phi", call, iterCnt);
//		write_debug(spF, "F", call, iterCnt);
////////////// BEBUG (end) ///////////////

		// normPhi = sqrt(normPhi);
//		UG_LOG(" ######################### (end sweep) ###################################  \n");
//		UG_LOG ("DamageFunctionUpdater: normPhi: "  << normPhi << "\n");	
	}

	m_lastNumIters = iterCnt;

	if(iterCnt >= maxIter){
		// UG_THROW("DamageFunctionUpdater: no convergence after " << iterCnt << " iterations");
		UG_LOG("DamageFunctionUpdater: no convergence after " << iterCnt << " iterations");
		return false;
	}


	UG_LOG ("DamageFunctionUpdater: maxPhi: "  << maxPhi << " after " <<iterCnt << " iterations\n");	

	return true;
}


template <typename TDomain>
void DamageFunctionUpdater<TDomain>::
set_debug(SmartPtr<GridFunctionDebugWriter<TDomain, CPUAlgebra> > spDebugWriter)
{
	m_spDebugWriter = spDebugWriter;
}

template <typename TDomain>
void DamageFunctionUpdater<TDomain>::
write_debug(SmartPtr<GridFunction<TDomain, CPUAlgebra> > spGF, std::string name, 
			int call, int iter)
{
	if(m_spDebugWriter.invalid()) return;

//	build name
	GridLevel gl = spGF->grid_level();
	std::stringstream ss;
	ss << "InDamageUpdate_" << name ;
	ss << "_call" << std::setfill('0') << std::setw(3) << call;
	if (iter >= 0) ss << "_iter" << std::setfill('0') << std::setw(3) << iter;
	ss << ".vec";

//	write
	m_spDebugWriter->set_grid_level(gl);
	m_spDebugWriter->write_vector(*spGF, ss.str().c_str());
}

template <typename TDomain>
void DamageFunctionUpdater<TDomain>::
write_stencil_matrix_debug(
			SmartPtr<GridFunction<TDomain, CPUAlgebra> > spGF, 
			std::string name, int call)
{
	if(m_spDebugWriter.invalid()) return;

//	build name
	GridLevel gl = spGF->grid_level();
	int numDoFs = spGF->num_dofs();
	std::stringstream ss;
	ss << "InDamageUpdate_" << name ;
	ss << "_call" << std::setfill('0') << std::setw(3) << call;
	//if (iter >= 0) ss << "_iter" << std::setfill('0') << std::setw(3) << iter;
	ss << ".mat";

	typedef CPUAlgebra::matrix_type TMat;			

	TMat A;
	A.resize_and_clear(numDoFs,numDoFs);

	for(size_t i = 0; i < m_vStencil.size(); ++i){
		for(size_t k = 0; k < m_vStencil[i].size(); ++k){
			const int j = m_vIndex[i][k];

			A(i,j) = m_vStencil[i][k];
		}
	}

//	write
	m_spDebugWriter->set_grid_level(gl);
	m_spDebugWriter->write_matrix(A, ss.str().c_str());

}



////////////////////////////////////////////////////////////////////////////////
// RelativeDensityUpdater
////////////////////////////////////////////////////////////////////////////////



template <typename TDomain>
std::vector<number> RelativeDensityUpdater<TDomain>::
solve(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spChi,
		SmartPtr<GridFunction<TDomain, CPUAlgebra> > spDrivingForce,
		const number betaStar, const number etaChiStar, 
		const number chiMin, const number dt, const int p,
		const number rho_target, const number MassTol)
{
	PROFILE_BEGIN_GROUP(RelativeDensityUpdater_solve, "Small Strain Mech");
	static int call = 0; call++;	



	////////////////////////////////////////////////////////////////////////////
	// rebuild after grid adaption
	////////////////////////////////////////////////////////////////////////////

	// get approximation space
	ConstSmartPtr<ApproximationSpace<TDomain> > approxSpace = spChi->approx_space();
	if(approxSpace != spDrivingForce->approx_space())
		UG_THROW("RelativeDensityUpdater<TDomain>::solve: expected same ApproximationSpace for f and psi0");

	// check revision counter if grid / approx space has changed since last call
	if(m_ApproxSpaceRevision != approxSpace->revision())
	{
		m_spElemSize = spChi->clone();
		m_spLaplaceChi = spChi->clone();
		m_spChiTrial = spChi->clone();

		// (re-)initialize setting
		InitLaplacianByPartialIntegration(m_spElemSize, m_vStencil, m_vIndex, m_quadRuleType, true);

		//	remember revision counter of approx space
		m_ApproxSpaceRevision = approxSpace->revision();

		// debug output
		write_stencil_matrix_debug(spChi, "STENCIL", call);
	}


	////////////////////////////////////////////////////////////////////////////
	// loop until maxiter
	////////////////////////////////////////////////////////////////////////////

	// TODO: better handling of h^2 in 3d
	number h2 = 0.0;
	for(size_t i = 0; i < m_vIndex.size(); ++i){
		const number vol = (*m_spElemSize)[i];
		h2 = std::max(h2, vol);
	}

	const int n = std::floor(6*betaStar / (etaChiStar * h2)) + 1;

	const number dt_chi = dt / n;

	int numBisect = 0; 

	// loop over inner (smaller) timesteps
	for(int j = 0; j < n; ++j)
	{
		////////////////////////////////////////////////////////////////////////////
		// compute parameter beta, eta
		////////////////////////////////////////////////////////////////////////////

		number min_p_chi = std::numeric_limits<number>::max();
		number max_p_chi = std::numeric_limits<number>::min();

		number int_g_p = 0.0, int_g = 0.0;
		for(size_t i = 0; i < m_vIndex.size(); ++i)
		{
			const number chi = (*spChi)[i];
			const number vol = (*m_spElemSize)[i];
			const number p_chi = (*spDrivingForce)[i];
			(*m_spLaplaceChi)[i] = Lambda(i, spChi);

			const number g = (chi - chiMin) * (1.0 - chi);

			int_g += vol * g;
			int_g_p += vol * g * p_chi;

			min_p_chi = std::min(min_p_chi, p_chi);
			max_p_chi = std::max(max_p_chi, p_chi);
		}

		const number p_w = int_g_p / int_g;

		const number beta = betaStar * p_w;
		const number eta_chi = etaChiStar * p_w;


		////////////////////////////////////////////////////////////////////////////
		// compute lagrange multiplier lambda
		////////////////////////////////////////////////////////////////////////////
/*
		number int_1 = 0.0, int_partPsi = 0.0, int_LaplaceChi = 0.0;
		for(size_t i = 0; i < m_vIndex.size(); ++i)
		{
			const number vol = (*m_spElemSize)[i];
			const number p_chi = (*spDrivingForce)[i];
			const number laplaceChi = (*m_spLaplaceChi)[i] = Lambda(i, spChi);

			int_1 += vol;
			int_partPsi += p_chi * vol;
			int_LaplaceChi += laplaceChi * vol;
		}

		const number lambda =  (1.0 / int_1) * (int_partPsi + beta * int_LaplaceChi);
*/
		////////////////////////////////////////////////////////////////////////////
		// update relative density (and theredy also driving force)
		////////////////////////////////////////////////////////////////////////////


		///////////// BEBUG (begin) ///////////////
//		write_debug(m_spLaplaceChi, "LaplaceChi", call, -1);
		///////////// BEBUG (end)  ///////////////
  	    number Llow = (min_p_chi - eta_chi);		// untere Grenze für Lagrangemultiplikator
		number Lhigh = (max_p_chi + eta_chi);	// obere Grenze für Lagrangemultiplikator
	 	number L_tr = 0.0; // why not: (Lhigh + Llow)/2.0; ???

	 	// bisection to find lagrange param
	 	while(true)
	 	{
	 		numBisect++;

			number Vol_Omega = 0.0;
	 		number rho_tr = 0.0;

	 		// Evolution equation for each element (trial step)
			for(size_t i = 0; i < m_vIndex.size(); ++i)
			{
				const number chi = (*spChi)[i];
				const number p_chi = (*spDrivingForce)[i];
				const number laplaceChi = (*m_spLaplaceChi)[i];
				const number vol = (*m_spElemSize)[i];

				number& chi_tr = (*m_spChiTrial)[i];

				// update relative density
				chi_tr = chi + (dt_chi / eta_chi) * (p_chi - L_tr + beta * laplaceChi);

				///////////// BEBUG (begin) ///////////////
				if(std::isfinite(chi_tr) == false){

					UG_LOG(" ###############  \n");
					UG_LOG("chi_tr     : " << chi_tr << "\n");
					UG_LOG("dt_chi     : " << dt_chi << "\n");
					UG_LOG("eta_chi    : " << eta_chi << "\n");
					UG_LOG("p_chi      : " << p_chi << "\n");
					UG_LOG("lambda     : " << L_tr << "\n");
					UG_LOG("beta       : " << beta << "\n");
					UG_LOG("laplaceChi : " << laplaceChi << "\n");
					UG_LOG("p_w        : " << p_w << "\n");
					UG_LOG(" ###############  \n");
					UG_THROW("Value for chi not finite, but: " << chi_tr);
				}
				////////////// BEBUG (end) ///////////////

				if(chi_tr > 1.0) chi_tr = 1.0;
				if(chi_tr < chiMin) chi_tr = chiMin;

				rho_tr += chi_tr * vol;
				Vol_Omega += vol;
			}

			rho_tr /= Vol_Omega;
			
			// Volume constraint
			// trial Volume
			if(fabs(rho_tr - rho_target) < MassTol) break; // V_tr = V_target

	        // Lagrange update		
			if(rho_tr > rho_target) Llow = L_tr;
			if(rho_tr < rho_target)	Lhigh = L_tr;
			
			L_tr = (Lhigh + Llow)/2.0;
		}


		// accept update
		for(size_t i = 0; i < m_vIndex.size(); ++i)
		{
			number& chi = (*spChi)[i];
			number& p_chi = (*spDrivingForce)[i];
			const number chi_tr = (*m_spChiTrial)[i];


			p_chi *= std::pow( (chi_tr / chi), p-1);

			chi = chi_tr;
		}
	}


	std::vector<number> vRes;

	vRes.push_back((number)n);
	vRes.push_back((number)numBisect);
	vRes.push_back((number)numBisect/(number)n);

	return vRes;
}




template <typename TDomain>
void RelativeDensityUpdater<TDomain>::
set_debug(SmartPtr<GridFunctionDebugWriter<TDomain, CPUAlgebra> > spDebugWriter)
{
	m_spDebugWriter = spDebugWriter;
}

template <typename TDomain>
void RelativeDensityUpdater<TDomain>::
write_debug(SmartPtr<GridFunction<TDomain, CPUAlgebra> > spGF, std::string name, 
			int call, int iter)
{
	if(m_spDebugWriter.invalid()) return;

//	build name
	GridLevel gl = spGF->grid_level();
	std::stringstream ss;
	ss << "InDensityUpdate_" << name ;
	ss << "_call" << std::setfill('0') << std::setw(3) << call;
	if (iter >= 0) ss << "_iter" << std::setfill('0') << std::setw(3) << iter;
	ss << ".vec";

//	write
	m_spDebugWriter->set_grid_level(gl);
	m_spDebugWriter->write_vector(*spGF, ss.str().c_str());
}

template <typename TDomain>
void RelativeDensityUpdater<TDomain>::
write_stencil_matrix_debug(
			SmartPtr<GridFunction<TDomain, CPUAlgebra> > spGF, 
			std::string name, int call)
{
	if(m_spDebugWriter.invalid()) return;

//	build name
	GridLevel gl = spGF->grid_level();
	int numDoFs = spGF->num_dofs();
	std::stringstream ss;
	ss << "InDensityUpdate_" << name ;
	ss << "_call" << std::setfill('0') << std::setw(3) << call;
	//if (iter >= 0) ss << "_iter" << std::setfill('0') << std::setw(3) << iter;
	ss << ".mat";

	typedef CPUAlgebra::matrix_type TMat;			

	TMat A;
	A.resize_and_clear(numDoFs,numDoFs);

	for(size_t i = 0; i < m_vStencil.size(); ++i){
		for(size_t k = 0; k < m_vStencil[i].size(); ++k){
			const int j = m_vIndex[i][k];

			A(i,j) = m_vStencil[i][k];
		}
	}

//	write
	m_spDebugWriter->set_grid_level(gl);
	m_spDebugWriter->write_matrix(A, ss.str().c_str());

}




////////////////////////////////////////////////////////////////////////////////
// Marking
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void MarkForAdaption_ValueRangeIndicator(
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spChi,
					IRefiner& refiner,
					number lowerValueToCoarsen, 
					number minValueToRefine, number maxValueToRefine,
					number upperValueToCoarsen,
					number maxJumpDiffToCoarsen,
					number minJumpDiffToRefine,
					int maxLevel)
{
	PROFILE_FUNC_GROUP("Small Strain Mech");

	static const int dim = TDomain::dim;
	typedef typename grid_dim_traits<dim>::element_type TElem; 
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;
	typedef typename TDomain::position_accessor_type	position_accessor_type;
	//position_accessor_type& aaPos = spChi->domain()->position_accessor();

	const int fct = 0; // \todo: generalize

	///////////////////////////
	// checks
	///////////////////////////

	if( minValueToRefine > maxValueToRefine ) UG_THROW("minValueToRefine > maxValueToRefine")
	if( lowerValueToCoarsen > minValueToRefine ) UG_THROW("lowerValueToCoarsen > minValueToRefine")
	if( upperValueToCoarsen < maxValueToRefine ) UG_THROW("upperValueToCoarsen < maxValueToRefine")



	///////////////////////////
	// marking
	///////////////////////////

	//	reset counter
	int numMarkedRefine = 0, numMarkedCoarse = 0;


	std::vector<MathVector<dim> > vCornerCoords;
	MathVector<dim> ElemCenter;

	std::vector<TElem*> vNeighbors;

	//	loop elements for marking
	const_iterator iter = spChi->template begin<TElem>();
	const_iterator iterEnd = spChi->template end<TElem>();
	for(; iter != iterEnd; ++iter)
	{
		//	get element
		TElem* elem = *iter;

		std::vector<DoFIndex> ind;
		if(spChi->inner_dof_indices(elem, fct, ind) != 1) UG_THROW("Wrong number dofs");
		const number val = DoFRef(*spChi, ind[0]);


		CollectSurfaceNeighbors(spChi, elem, vNeighbors);

	//	check whether all children of e of type TElem have similar values
		bool neighborsHaveLargeJump = false;
		for(size_t i = 0; i < vNeighbors.size(); ++i){

			TElem* neighbor = vNeighbors[i];

			std::vector<DoFIndex> ind;
			if(spChi->inner_dof_indices(neighbor, fct, ind) != 1) UG_THROW("Wrong number dofs");
			const number valNeighbor = DoFRef(*spChi, ind[0]);

			if(fabs(val - valNeighbor) > minJumpDiffToRefine)
			{
				neighborsHaveLargeJump = true;
				break;
			}
		}

		//	marks for refinement
		if(      (val >= minValueToRefine && val <= maxValueToRefine) 
		     ||  (neighborsHaveLargeJump) )
		{

			if(spChi->dd()->multi_grid()->get_level(elem) < maxLevel)
			{
				refiner.mark(elem, RM_REFINE);
				numMarkedRefine++;
			}

			continue;
		}

		//	marks for coarsening
		if( val < lowerValueToCoarsen || val > upperValueToCoarsen)
		{

		//	get the parent
			TElem* parent = dynamic_cast<TElem*>(spChi->dd()->multi_grid()->get_parent(elem));
			if(parent){

			//	check whether all children of e of type TElem have similar values
				bool allChildHaveSimilarValue = true;
				size_t numChildren = spChi->dd()->multi_grid()->template num_children<TElem>(parent);
				for(size_t i = 0; i < numChildren; ++i){

					TElem* child = spChi->dd()->multi_grid()->template get_child<TElem>(parent, i);

					// check if covered: cannot coarsen covered children (i.e. with children themselves)
					if(spChi->dd()->multi_grid()->template num_children<TElem>(child) > 0){
						allChildHaveSimilarValue = false;
						break;						
					}

					// check if values are the same for all (leaf-)children
					std::vector<DoFIndex> ind;
					if(spChi->inner_dof_indices(child, fct, ind) != 1) UG_THROW("Wrong number dofs");
					const number valChild = DoFRef(*spChi, ind[0]);
					if(fabs(val - valChild) > maxJumpDiffToCoarsen)
					{
						allChildHaveSimilarValue = false;
						break;
					}
				}
			

				if(allChildHaveSimilarValue){					
					refiner.mark(elem, RM_COARSEN);
					numMarkedCoarse++;
				}
			}

		}

	}

	///////////////////////////
	// print infos
	///////////////////////////
	UG_LOG("MarkForAdaption_ValueRangeIndicator:\n")
#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1)
	{
		UG_LOG("  >>> Marked on Proc "<<pcl::ProcRank()<<": refine: " << numMarkedRefine << ", coarsen : " << numMarkedCoarse << "\n");
		pcl::ProcessCommunicator com;
		int numMarkedRefineLocal = numMarkedRefine, numMarkedCoarseLocal = numMarkedCoarse;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
		numMarkedCoarse = com.allreduce(numMarkedCoarseLocal, PCL_RO_SUM);
	}
#endif

	UG_LOG("  >>> Marked refine: " << numMarkedRefine << ", coarsen: " << numMarkedCoarse << "\n" );
}


template<typename TDomain>
void MarkDamage(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0,
					IRefiner& refiner,
					number minValueToRefine, number maxValueToCoarsen,
					int maxLevel, 
					const std::vector<MathVector<TDomain::dim,number>* >& vCenter, 
					const std::vector<number>& vRadius)
{
	PROFILE_FUNC_GROUP("Small Strain Mech");

	static const int dim = TDomain::dim;
	typedef typename grid_dim_traits<dim>::element_type TElem; 
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;
	typedef typename TDomain::position_accessor_type	position_accessor_type;
	position_accessor_type& aaPos = spF->domain()->position_accessor();

	const int fct = 0; // \todo: generalize

	///////////////////////////
	// checks
	///////////////////////////

	if( (minValueToRefine < maxValueToCoarsen) )
		UG_THROW("Coarsen threshold must be smaller than refine threshold")

	if( !(vRadius.size() == vCenter.size()) )
		UG_THROW("Num radius and center must match")
	

	///////////////////////////
	// marking
	///////////////////////////

	//	we'll compare against the square radius.
	std::vector<number> vRadiusSq;
	for(size_t i = 0; i < vRadius.size(); ++i)
		vRadiusSq.push_back( vRadius[i] * vRadius[i]);

	std::vector<MathVector<dim> > vCircleCenter;
	for(size_t i = 0; i < vCenter.size(); ++i)
		vCircleCenter.push_back( *(vCenter[i]));


	//	reset counter
	int numMarkedRefine = 0, numMarkedCoarse = 0;


	std::vector<MathVector<dim> > vCornerCoords;
	MathVector<dim> ElemCenter;

	//	loop elements for marking
	const_iterator iter = spF->template begin<TElem>();
	const_iterator iterEnd = spF->template end<TElem>();
	for(; iter != iterEnd; ++iter)
	{
		//	get element
		TElem* elem = *iter;

		std::vector<DoFIndex> ind;
		if(spF->inner_dof_indices(elem, fct, ind) != 1) UG_THROW("Wrong number dofs");
		const number val = DoFRef(*spF, ind[0]) * DoFRef(*spPsi0, ind[0]);

		// check for extra refinement
		bool bInCircle = false;
		if(!vRadiusSq.empty()){
			CollectCornerCoordinates(vCornerCoords, *elem, aaPos);
			AveragePositions<dim>(ElemCenter, vCornerCoords);

			for(size_t i = 0; i < vRadiusSq.size(); ++i){
				if(VecDistanceSq(vCircleCenter[i], ElemCenter) <= vRadiusSq[i]){
					bInCircle = true;
				}
			}
		}

		//	marks for refinement
		if( val > minValueToRefine || bInCircle)
		{
			if(spF->dd()->multi_grid()->get_level(elem) < maxLevel)
			{
				refiner.mark(elem, RM_REFINE);
				numMarkedRefine++;
			}

			continue;
		}

		//	marks for coarsening
		if( val < maxValueToCoarsen)
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
		UG_LOG("  >>> Marked on Proc "<<pcl::ProcRank()<<": refine: " << numMarkedRefine << ", coarsen : " << numMarkedCoarse << "\n");
		pcl::ProcessCommunicator com;
		int numMarkedRefineLocal = numMarkedRefine, numMarkedCoarseLocal = numMarkedCoarse;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
		numMarkedCoarse = com.allreduce(numMarkedCoarseLocal, PCL_RO_SUM);
	}
#endif

	UG_LOG("  >>> Marked refine: " << numMarkedRefine << ", coarsen: " << numMarkedCoarse << "\n" );
}



template<typename TDomain>
void MarkDamage_OLD_AND_DEPRECATED(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0,
					IRefiner& refiner,
					number refineFrac, number coarseFrac, 
					number avgRefineFactor, number avgCoarsenFactor,
					int maxLevel)
{
	PROFILE_FUNC_GROUP("Small Strain Mech");

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
std::vector<number> DamageStatistic(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
										SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0)
{
	PROFILE_FUNC_GROUP("Small Strain Mech");

	static const int dim = TDomain::dim;
	typedef typename grid_dim_traits<dim>::element_type TElem; 
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;
	const int fct = 0;

	///////////////////////////
	// statistic
	///////////////////////////

	number max_fPsi0 = std::numeric_limits<number>::min();
	number min_fPsi0 = std::numeric_limits<number>::max();
	number sum_fPsi0 = 0.0;

	number max_l2_fPsi0 = std::numeric_limits<number>::min();
	number min_l2_fPsi0 = std::numeric_limits<number>::max();
	number sum_l2_fPsi0 = 0.0;

	number max_l_fPsi0 = std::numeric_limits<number>::min();
	number min_l_fPsi0 = std::numeric_limits<number>::max();
	number sum_l_fPsi0 = 0.0;

	// loop all elems 
	const_iterator iter = spF->template begin<TElem>();
	const_iterator iterEnd = spF->template end<TElem>();

	//	loop elements for marking
	size_t numElem = 0;
	for(; iter != iterEnd; ++iter)
	{
		//	get element
		TElem* elem = *iter; ++numElem;

		std::vector<DoFIndex> ind;
		if(spF->inner_dof_indices(elem, fct, ind) != 1) UG_THROW("Wrong number dofs");

		const number f =  DoFRef(*spF, ind[0]);
		const number psi0 = DoFRef(*spPsi0, ind[0]);

		max_fPsi0 = std::max(max_fPsi0, f*psi0);
		min_fPsi0 = std::min(min_fPsi0, f*psi0);
		sum_fPsi0 += f*psi0;

		const number l2 = ElementDiameterSq(*elem, *(spF->domain()) );

		max_l2_fPsi0 = std::max(max_l2_fPsi0, l2*f*psi0);
		min_l2_fPsi0 = std::min(min_l2_fPsi0, l2*f*psi0);
		sum_l2_fPsi0 += l2*f*psi0;

		const number l = std::sqrt(l2);

		max_l_fPsi0 = std::max(max_l_fPsi0, l*f*psi0);
		min_l_fPsi0 = std::min(min_l_fPsi0, l*f*psi0);
		sum_l_fPsi0 += l*f*psi0;
	}
	number avg_fPsi0 = sum_fPsi0 / numElem;
	number avg_l2_fPsi0 = sum_l2_fPsi0 / numElem;
	number avg_l_fPsi0 = sum_l_fPsi0 / numElem;


	std::vector<number> vRes;
	vRes.push_back(max_fPsi0);
	vRes.push_back(min_fPsi0);
	vRes.push_back(avg_fPsi0);

	vRes.push_back(max_l2_fPsi0);
	vRes.push_back(min_l2_fPsi0);
	vRes.push_back(avg_l2_fPsi0);

	vRes.push_back(max_l_fPsi0);
	vRes.push_back(min_l_fPsi0);
	vRes.push_back(avg_l_fPsi0);

	return vRes;
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




template<typename TDomain>
std::vector<number> MinMaxElementDiameter(SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF)
{
	PROFILE_FUNC_GROUP("Small Strain Mech");

	static const int dim = TDomain::dim;
	typedef typename grid_dim_traits<dim>::element_type TElem; 
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

	typename TDomain::grid_type& grid = *(spF->domain()->grid());
	typename TDomain::position_accessor_type& aaPos = spF->domain()->position_accessor();

	// loop all elems 
	const_iterator iter = spF->template begin<TElem>();
	const_iterator iterEnd = spF->template end<TElem>();

	//	loop elements for marking
	number max = 0.0;
	number min = std::numeric_limits<number>::max();
	for(; iter != iterEnd; ++iter){
		const number size = ElementDiameterSq(grid, aaPos, *iter);
		max = std::max(max, size);
		min = std::min(min, size);
	}

#ifdef UG_PARALLEL
	// share value between all procs
	pcl::ProcessCommunicator com;
	max = com.allreduce(max, PCL_RO_MAX);
	min = com.allreduce(min, PCL_RO_MIN);
#endif

	std::vector<number> vRes;

	vRes.push_back( std::sqrt(min) );
	vRes.push_back( std::sqrt(max) );

	return vRes;

}




} // end namespace SmallStrainMechanics
}// namespace ug

#endif /* __SMALL_STRAIN_MECHANICS__DAMAGE_IMPL_H_ */
