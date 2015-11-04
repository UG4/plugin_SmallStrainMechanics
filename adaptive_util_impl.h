
#ifndef ADAPTIVE_UTIL_IMPL_H_
#define ADAPTIVE_UTIL_IMPL_H_

#include "lib_grid/algorithms/refinement/refiner_interface.h"
#include "adaptive_util.h"

namespace ug{
namespace SmallStrainMechanics{

template<typename TDomain, typename TElem>
void plastic_ip_elem(bool& bPlasticIPs, TElem* elem,
		const LocalVector& locU, SmartPtr<TDomain> dom,
		SmallStrainMechanicsElemDisc<TDomain>& elemDisc)
{
	static const int dim = TDomain::dim;

	//	get vertices and extract corner coordinates
	typedef typename TDomain::position_accessor_type position_accessor_type;
	position_accessor_type& aaPos = dom->position_accessor();
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	const size_t numVertices = elem->num_vertices();
	for (size_t i = 0; i < numVertices; ++i) {
		coCoord[i] = aaPos[elem->vertex(i)];
	};

	int quadOrder = elemDisc.get_quad_order();

	//	update the geometry for the Element elem
	DimFEGeometry<dim> geo;
   	try{
		geo.update(elem, &(coCoord[0]),
				LFEID(LFEID::LAGRANGE, dim, 1), quadOrder);
	}
   	UG_CATCH_THROW("plastic_ip_elem:"
		" Cannot update Finite Element Geometry.");

   	//	get the material law
   	if (elemDisc.get_material_law().invalid())
   		UG_THROW("No material law set in plastic_ip_elem \n");

	SmartPtr<IMaterialLaw<TDomain> > spMatLaw = elemDisc.get_material_law();

	//  pointer to internal variable of current elem
	spMatLaw->internal_vars(elem);

	MathMatrix<dim, dim> GradU;
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//	get displacementGradient (GradU)
		spMatLaw->template DisplacementGradient<DimFEGeometry<dim> >(GradU, ip, geo, locU);

		number gamma = spMatLaw->plastic_multiplier(ip, GradU);
		if (gamma > 0.0)
		{
			bPlasticIPs = true; return;
		}
		if (gamma == 0.0)
		{
			bPlasticIPs = false;
		}
		if (gamma < 0.0)
			UG_THROW("gamma: " << gamma << "in plastIP_elem \n");
	}
}

template <typename TDomain, typename TAlgebra>
void MarkForAdaption_PlasticElem(IRefiner& refiner,
     GridFunction<TDomain, TAlgebra>& u,
     SmallStrainMechanicsElemDisc<TDomain>& elemDisc)
{
//	types
	typedef GridFunction<TDomain, TAlgebra> TFunction;
	typedef typename TFunction::element_type element_type;
	typedef typename DoFDistribution::traits<element_type>::const_iterator const_iterator;

	ConstSmartPtr<DoFDistribution> dd = u.dof_distribution();

	int numMarkedRefine = 0;
	bool bPlasticIPs = false;
	const_iterator iter = dd->template begin<element_type>();
	const_iterator iterEnd = dd->template end<element_type>();

	// 	local indices and local algebra
	LocalIndices ind; LocalVector locU;

//	loop elements for marking
	for(; iter != iterEnd; ++iter)
	{
		element_type* elem = *iter;

		// 	get global indices
		u.indices(elem, ind);

		// 	adapt local algebra
		locU.resize(ind);

		//	local vector extract -> locU
		GetLocalVector(locU, u);

		//	reset boolean
		bPlasticIPs = false;

		//	check, if plastic yielding occurs
		//	at one integration point in the element
		plastic_ip_elem<TDomain, element_type> (bPlasticIPs, elem, locU,
				u.domain(), elemDisc);

		//	marks for refinement
		if(bPlasticIPs)
		{
			refiner.mark(elem, RM_REFINE);
			numMarkedRefine++;
		}
	}

#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1)
	{
		UG_LOG("  +++ Marked for refinement on Proc "<<pcl::ProcRank()<<": " << numMarkedRefine << " Elements.\n");
		pcl::ProcessCommunicator com;
		int numMarkedRefineLocal = numMarkedRefine;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
	}
#endif

UG_LOG("  +++ Marked for refinement: " << numMarkedRefine << " Elements.\n");

}



} //end of namespace SmallStrainMechanics
} //end of namespace ug

#endif /* ADAPTIVE_UTIL_IMPL_H_ */
