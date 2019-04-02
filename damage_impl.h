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

namespace ug{
namespace SmallStrainMechanics{


template <typename TDomain>
void DamageFunctionUpdater<TDomain>::
AveragePositions(	MathVector<dim>& vCenter, 
					const std::vector<MathVector<dim> >& vCornerCoords)
{
	vCenter = vCornerCoords[0];
	for(size_t j = 1; j < vCornerCoords.size(); ++j)
	{
		vCenter += vCornerCoords[j];
	}
	vCenter *= 1./(number)( vCornerCoords.size());
}


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
	AveragePositions(vDistance.back(), vCornerCoords);
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

/*				UG_LOG(" ++++++ Distances  Elem "  << i << "\n")
		for (int k = 0; k < numNeighbors; ++k)
			UG_LOG(" ++"  << k << ": "  << vDistance[k] <<"\n")
		UG_LOG(" +++ DLambda "  << m_vDLambda[i] << "\n")
*/
		
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


template <typename TDomain>
void DamageFunctionUpdater<TDomain>::
init_ByPartIntegral(	
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
					SmartPtr<GridFunction<TDomain, CPUAlgebra> > spPsi0)
{
	PROFILE_BEGIN_GROUP(DamageFunctionUpdater_init_ByPartIntegral, "Small Strain Mech");

	const size_t fct = 0;

	// weights (for Simpson's rule)
	const number sideWeight = 4.0/6.0;
	const number vertexWeight = 1.0/6.0;

//			const number sideWeight = 1.0;
//			const number vertexWeight = 0.0;


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
	m_vIndex.resize(numIndex);
	m_vStencil.resize(numIndex);

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
		m_vIndex[i].clear();
		m_vIndex[i].push_back(i);				

		// reset stencil
		m_vStencil[i].clear(); 
		m_vStencil[i].resize(1, 0.0);

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

		// element midpoint
		AveragePositions(ElemMidPoint, vCornerCoords);

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
						AveragePositions(Distance, vNbrCornerCoords);
						VecScaleAdd(Distance, 1.0, ElemMidPoint, -1.0, Distance);

						const number cpl = (sideWeight) * VecDot(Normal, Distance) / ( VecTwoNormSq(Distance) * vol);
						m_vStencil[i].push_back( -cpl );
						m_vStencil[i][0]   +=  cpl;

						std::vector<DoFIndex> ind;
						if(spF->inner_dof_indices(neighborElem, fct, ind) != 1) UG_THROW("Wrong number dofs");
						m_vIndex[i].push_back(ind[0][0]);

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
					m_vStencil[i].push_back( -cpl );
					m_vStencil[i][0]   +=  cpl;

					std::vector<DoFIndex> ind;
					if(spF->inner_dof_indices(elem, fct, ind) != 1) UG_THROW("Wrong number dofs");
					m_vIndex[i].push_back(ind[0][0]);
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
					AveragePositions(Distance, vNbrCornerCoords);
					VecScaleAdd(Distance, 1.0, ElemMidPoint, -1.0, Distance);

					const number cpl = (sideWeight / numConstrained) * VecDot(Normal, Distance) / ( VecTwoNormSq(Distance) * vol);
					m_vStencil[i].push_back( -cpl );
					m_vStencil[i][0]   +=  cpl;

					std::vector<DoFIndex> ind;
					if(spF->inner_dof_indices(neighborElem, fct, ind) != 1) UG_THROW("Wrong number dofs");
					m_vIndex[i].push_back(ind[0][0]);
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
					AveragePositions(Distance, vNbrCornerCoords);
					VecScaleAdd(Distance, 1.0, ElemMidPoint, -1.0, Distance);

					const number cpl = sideWeight * VecDot(Normal, Distance) / ( VecTwoNormSq(Distance) * vol);
					m_vStencil[i].push_back( -cpl );
					m_vStencil[i][0]   +=  cpl;

					std::vector<DoFIndex> ind;
					if(spF->inner_dof_indices(neighborElem, fct, ind) != 1) UG_THROW("Wrong number dofs");
					m_vIndex[i].push_back(ind[0][0]);
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
						AveragePositions(Distance, vNbrCornerCoords);
						VecScaleAdd(Distance, 1.0, ElemMidPoint, -1.0, Distance);


						//UG_LOG("Elem: "<<i<<", side: "<<s<<", vrt: "<<vrt<<", Distance: "<<Distance<<", Normal: "<<Normal<<"\n");

						if(dim == 3)
							UG_THROW("This implementation is 2d only, currently. Handle vertex neighbors properly in 3d...");

						std::vector<DoFIndex> ind;
						if(spF->inner_dof_indices(neighborElem, fct, ind) != 1) UG_THROW("Wrong number dofs");
						const size_t j = ind[0][0];

						size_t k = 0;
						for(; k < m_vIndex[i].size(); ++k){
							if(m_vIndex[i][k] == j) 
								break;
						}
						if(k == m_vIndex[i].size()) 
						{
							m_vIndex[i].push_back(j);
							m_vStencil[i].push_back(0.0);
							// k == m_vIndex[i].size();
						}

						const number cpl = (vertexWeight / numDiagElem) * VecDot(Normal, Distance) / ( VecTwoNormSq(Distance) * vol);
						m_vStencil[i][k]   -= cpl;
						m_vStencil[i][0]   +=  cpl;
					}
				}
			} // end side loop

		} // end diag-elem block

		//	end marking
		grid.end_marking();

	} // end element loop

	static int call = 0; call++;
	write_stencil_matrix_debug(spF, "Stencil", call);

}

template <typename TDomain>
number DamageFunctionUpdater<TDomain>::
Lambda_ByPartIntegral(size_t i, SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF) 
{

	number res = 0.0;
	for (size_t j = 0; j < m_vIndex[i].size(); ++j){

		res += m_vStencil[i][j] * (*spF)[ m_vIndex[i][j] ];
	}

	return res;
}


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
		init_ByPartIntegral(spF, spPsi0);

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
		(*spLambdaOld)[i] = DLambda_ByPartIntegral(i);
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
			(*spLambdaOld)[i] = Lambda_ByPartIntegral(i, spF);

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
				f = f - dampNewton * (phi / (psi0 - beta * (lambda + f * DLambda_ByPartIntegral(i)) ));

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
					UG_LOG("DLambda: " << DLambda_ByPartIntegral(i) << "\n");
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


	UG_LOG ("DamageFunctionUpdater: normPhi: "  << normPhi << " after " <<iterCnt << " iterations\n");	

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



template<typename TDomain>
void MarkDamage(	SmartPtr<GridFunction<TDomain, CPUAlgebra> > spF,
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


} // end namespace SmallStrainMechanics
}// namespace ug

#endif /* __SMALL_STRAIN_MECHANICS__DAMAGE_IMPL_H_ */
