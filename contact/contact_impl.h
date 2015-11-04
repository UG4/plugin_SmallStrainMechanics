
#ifndef CONTACT_IMPL_H_
#define CONTACT_IMPL_H_

// other ug4 modules
#include "common/common.h"
#include "lib_disc/common/geometry_util.h"

#include "contact.h"

namespace ug{
namespace SmallStrainMechanics{

template <int dim> struct face_type_traits
{
    typedef void face_type0;
	typedef void face_type1;
	typedef void DimFEGeo;
};

template <> struct face_type_traits<1>
{
    typedef ReferenceVertex face_type0;
	typedef ReferenceVertex face_type1;
	typedef DimFEGeometry<1, 1> DimFEGeo;
};

template <> struct face_type_traits<2>
{
    typedef ReferenceEdge face_type0;
	typedef ReferenceEdge face_type1;
	typedef DimFEGeometry<2, 1> DimFEGeo;
};

template <> struct face_type_traits<3>
{
    typedef ReferenceTriangle face_type0;
	typedef ReferenceQuadrilateral face_type1;
	typedef DimFEGeometry<3, 2> DimFEGeo;
};

template <typename TDomain>
template <typename TElem>
void SmallStrainMechanicsElemDisc<TDomain>::contact_forces_elem_ips_avg(
		LocalVector& locForce, GridObject* side,
		TElem* elem, const MathVector<dim> sideCoPos[], int numElemCorners,
		const LocalVector& locU, vector<DoFIndex> vActiveSetLoc)
{
	typedef typename face_type_traits<dim>::face_type0 face_type0;
	typedef typename face_type_traits<dim>::face_type1 face_type1;
	typedef typename face_type_traits<dim>::DimFEGeo sideGeo;

	ReferenceObjectID sideRoid = side->reference_object_id();
	sideGeo geo(sideRoid, 3, LFEID(LFEID::LAGRANGE, dim, 1));

	//	prepare geometry for type and order
	try{
		geo.update_local(sideRoid, LFEID(LFEID::LAGRANGE, dim, 1), 3);
	}
	UG_CATCH_THROW("SmallStrainMechanicsElemDisc::contact_forces_elem_ips_avg:"
			" Cannot update local values of finite element geometry.");

	try{
		geo.update(side, sideCoPos, LFEID(LFEID::LAGRANGE, dim, 1), 3);
	}
	UG_CATCH_THROW("SmallStrainMechanicsElemDisc::contact_forces_elem_ips_avg:"
					" Cannot update finite element geometry.");

	MathVector<dim> normal, normalizedNormal;
	if (numElemCorners == dim)
		ElementNormal<face_type0, dim>(normal, sideCoPos);
	else
		ElementNormal<face_type1, dim>(normal, sideCoPos);

	number normOfNormal = VecLength(normal);
	VecScale(normalizedNormal, normal, 1.0/normOfNormal);

	MathMatrix<dim, dim> GradU, sigma;
	for(size_t ip = 0; ip < geo.num_ip(); ++ip) // loop ip
	{
		//	get displacementGradient (GradU), linearized strain tensor (eps)
		//	and Cauchy stress tensor sigma at ip

		m_spMatLaw->template DisplacementGradient<sideGeo>(GradU, ip, geo, locU);
		m_spMatLaw->stressTensor(sigma, ip, GradU);

		/*for (size_t i = 0; i < (size_t) dim; ++i){
			for (size_t J = 0; J < (size_t) dim; ++J)
			{
				GradU[i][J] = 0.0;
				for (size_t a = 0; a < geo.num_sh(); ++a) // loop shapes
					GradU[i][J] += geo.global_grad(ip, a)[J] * locU(i, a);
			}
		}

		for (size_t i = 0; i < (size_t) dim; ++i)
			for (size_t J = 0; J < (size_t) dim; ++J)
				eps[i][J] = 0.5 * (GradU[i][J] + GradU[J][i]);

		TensContract4(sigma, m_ElastTensorFunct, eps);*/

		number sigma_n = 0.0;
		for (size_t i = 0; i < (size_t) dim; ++i)
			for (size_t j = 0; j < (size_t) dim; ++j)
				sigma_n += sigma[i][j] * normalizedNormal[i] * normalizedNormal[j];

		//	compute normal component of sigma; sigma_n = sigma_ij * normal_i * normal_j
		//	'vActiveSetLoc' is a vector, which stores all active local (dof,fct)-pairs
		//	of the SideElem 'side'
		for (vector<DoFIndex>::iterator itActiveInd = vActiveSetLoc.begin();
				itActiveInd < vActiveSetLoc.end(); ++itActiveInd)
			locForce((*itActiveInd)[1], (*itActiveInd)[0]) += sigma_n * geo.weight(ip);

	} // end(ip)

	//	FOR 3D-FULL-HEXAEDER-ELEMENTS instead of sideElems it is possible
	//	to consider plasticity by means of plastic ip's via
	//	'm_pElemData = &m_aaElemData[static_cast<TElem*>(elem)];'
}

template <typename TDomain>
template <typename TSide, typename TElem>
void SmallStrainMechanicsElemDisc<TDomain>::contact_forces_elem_midpoint(
		LocalVector& locForce, TSide* side,
		TElem* elem, const MathVector<dim> sideCoPos[],
		const LocalVector& locU, vector<DoFIndex> vActiveSetLoc)
{
	//	TODO: inclusion of plastic variables, which are defined at the elements ip`s
	typedef typename face_type_traits<dim>::face_type0 face_type0;
	typedef typename face_type_traits<dim>::face_type1 face_type1;
	SmartPtr<TDomain> dom = this->domain();
	typename domain_type::position_accessor_type& aaPos = dom->position_accessor();

	//	reference object type
	ReferenceObjectID sideRoid = side->reference_object_id();

//	some storage
	MathMatrix<dim, dim-1> JTInv;
	vector<MathVector<dim-1> > vLocalGrad;
	vector<MathVector<dim> > vGlobalGrad;
	vector<MathVector<dim> > vCorner;

//	get trial space
	const LocalShapeFunctionSet<dim-1>& lsfs =
			LocalFiniteElementProvider::get<dim-1>(sideRoid, LFEID(LFEID::LAGRANGE, dim, 1));

//	create a reference mapping
	DimReferenceMapping<dim-1, dim>& refMap
		= ReferenceMappingProvider::get<dim-1, dim>(sideRoid);

//	number of shape functions
	const size_t numSH = lsfs.num_sh();
	vLocalGrad.resize(numSH);
	vGlobalGrad.resize(numSH);

//	get local Mid Point (for ROID_EDGE, ROID_QUADRILATERAL, ROID_TRIANGLE)
	if (dim-1 < 1 || dim-1 > 2)
		UG_THROW("SmallStrainMechanicsElemDisc::contact_forces_elem_midpoint:"
				"side dimension " << dim-1 << " is not supported!");

	MathVector<dim-1> localIP;
	for (int i = 0; i < dim-1; ++i)
		localIP[i] = 0.5;

	if (sideRoid == ROID_TRIANGLE){
		localIP[0] = 1.0/3.0; localIP[1] = 1.0/3.0;
	}

//	evaluate reference gradient at local midpoint
	lsfs.grads(&vLocalGrad[0], localIP);

//	get corners of element
	CollectCornerCoordinates(vCorner, *side, aaPos);

//	update mapping
	refMap.update(&vCorner[0]);

//	compute jacobian
	refMap.jacobian_transposed_inverse(JTInv, localIP);

//	call 'ElementNormal' depending on elem-type
//	faces have dim corners in 1d, 2d
//	in 3d they have dim corners (triangle) or dim+1 corners (quadrilateral)
	MathVector<dim> normal, normalizedNormal;

	int numElemCorners = (int)vCorner.size();
	if (numElemCorners == dim)
		ElementNormal<face_type0, dim>(normal, sideCoPos);
	else
		ElementNormal<face_type1, dim>(normal, sideCoPos);

	number normOfNormal = VecLength(normal);
	VecScale(normalizedNormal, normal, 1.0/normOfNormal);

//	compute gradient at mid point by summing contributions of all shape fct
	MathMatrix<dim, dim> GradU, eps, sigma;
	GradU = 0.0;

	sigma = 0.0;
	/*for(size_t sh = 0; sh < numSH; ++sh)
	{
	//	get global Gradient
		MatVecMult(vGlobalGrad[sh], JTInv, vLocalGrad[sh]);
		for(size_t fct = 0; fct < locU.num_fct(); ++fct)
			for(size_t i = 0; i < (size_t) dim; ++i)
				GradU[fct][i] += vGlobalGrad[sh][i] * locU(fct, sh);
	}

	//	get strain tensor 'eps' at local midPoint-IP
	for(size_t fct = 0; fct < locU.num_fct(); ++fct)
		for (size_t i = 0; i < (size_t) dim; ++i)
			eps[fct][i] = 0.5 * (GradU[fct][i] + GradU[i][fct]);

	//	get Cauchy stress tensor 'sigma'
	TensContract4(sigma, m_ElastTensorFunct, eps);*/

	//	compute normal component -sigma_n
	number sigma_n = 0.0;
	for (size_t i = 0; i < (size_t) dim; ++i)
		for (size_t j = 0; j < (size_t) dim; ++j)
			sigma_n += sigma[i][j] * normalizedNormal[i] * normalizedNormal[j];

/*	UG_LOG("sigma_n: " << sigma_n << "\n");
	MathVector<dim> sigma_n_normal;
	VecScale(sigma_n_normal, normalizedNormal, sigma_n);
	UG_LOG("sigma_n_normal: " << sigma_n_normal << "\n");*/

	//	add contributions to contactForces for active indices
	for (vector<DoFIndex>::iterator itAcvtiveInd = vActiveSetLoc.begin();
			itAcvtiveInd < vActiveSetLoc.end(); ++itAcvtiveInd)
		locForce((*itAcvtiveInd)[1], (*itAcvtiveInd)[0]) = sigma_n; //sigma_n_normal[activeLocFct];

}

template <typename TDomain, typename TGridFunction>
ContactSmallStrainMechanics<TDomain, TGridFunction>::ContactSmallStrainMechanics(
		SmartPtr<SmallStrainMechanicsElemDisc<TDomain> > spMechElemDisc)
{
	m_spMechElemDisc = spMechElemDisc;
}

template <typename TDomain, typename TGridFunction>
template <typename TSideElem, typename TIterator>
void ContactSmallStrainMechanics<TDomain, TGridFunction>::contact_forces_elem(
		TIterator iterBegin,
		TIterator iterEnd,
		const TGridFunction& u,
		TGridFunction& contactForce,
		vector<DoFIndex> vActiveSetGlob)
{
// 	check if at least an element exists, else return
	if(iterBegin == iterEnd) return;

	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::grid_type grid_type;
	typedef typename TGridFunction::element_type element_type;
	typename grid_type::template traits<element_type>::secure_container associatedElems;
	typename domain_type::position_accessor_type& aaPos
				= contactForce.domain()->position_accessor();

	grid_type& grid = *contactForce.domain()->grid();
	vector<DoFIndex> vActiveSetLoc;

	// 	local indices and local algebra
	LocalIndices ind; LocalVector locU, locContactForce;

//	loop over all elements on activeSubset si
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		TSideElem* sideElem = *iter;

		//	get all associated elements to 'elem' of activeSubset
		grid.associated_elements(associatedElems, sideElem);

		if (associatedElems.size() > 1)
			UG_THROW("Temporarily it is only possible to compute contact Forces "
					"for an element of dimension " << dim-1 << " with exactly "
					"one associated element of dimension " << dim << " in "
					"ContactSmallStrainMechanics::contactForcesElem \n");

		// 	get global indices
		u.indices(*iter, ind);

		// 	adapt local algebra
		locU.resize(ind); locContactForce.resize(ind);

		vActiveSetLoc.resize(0);

		//	create local activeSet-vector out of global activeSet-vector
		for (vector<DoFIndex>::iterator itActiveInd = vActiveSetGlob.begin();
				itActiveInd < vActiveSetGlob.end(); ++itActiveInd)
		{
			for(size_t fct=0; fct < locU.num_all_fct(); ++fct){
				for(size_t dof=0; dof < locU.num_all_dof(fct); ++dof)
				{
					size_t globIndex = ind.index(fct, dof);
					size_t globComp = ind.comp(fct, dof);

					//	check if local (fct,dof) corresponds to a global activeSetMultiIndex
					if ((*itActiveInd)[0] == globIndex &&
							(*itActiveInd)[1] == globComp)
					{
						UG_LOG("active multiIndex (ind,comp) in "
								"ContactSmallStrainMechanics::contactForces: ("
								<< globIndex << "," << globComp << ") \n");

						//	store local (dof,fct)-pairs which are active
						//	in the current element
						DoFIndex activeMultiIndexLoc(dof, fct);

						UG_LOG("active locU (dof,fct) in "
								"ContactSmallStrainMechanics::contactForces: ("
								<< dof << "," << fct << ") \n");

						//	create list of active MultiIndex-pairs (local)
						vActiveSetLoc.push_back(activeMultiIndexLoc);
					}
				}
			}
		} // end (vActiveSet-loop)

		if (vActiveSetLoc.size() != 0.0)
		{
			UG_LOG("activeDoFs in Elem" << *iter << "in "
					"ContactSmallStrainMechanics::contactForces: "
					<< vActiveSetLoc.size() << "\n");

			//	reset contribution of this element
			locContactForce = 0.0;

			//	local vector extract -> locU
			GetLocalVector(locU, u);

			//	storage for corner coordinates
			vector<MathVector<dim> > vCorner;
			MathVector<dim> sideCoPos[dim+1];

			//	reference object type
			ReferenceObjectID roid = sideElem->reference_object_id();

			const DimReferenceElement<dim-1>& rRefElem
					= ReferenceElementProvider::get<dim-1>(roid);

			//	get corner coordinates
			CollectCornerCoordinates(vCorner, *sideElem, aaPos);

			//	here the ordering of the corners in the reference element is exploited
			//	in order to compute the outer normal later on
			int numElemCorners = (int)vCorner.size();

			for (int i = 0; i < numElemCorners; ++i)
				sideCoPos[i] = vCorner[rRefElem.id(dim-1, 0, 0, i)];

			//	pass the local indices (dof,fct) which are active
			//	to elementwise computation of contactForces
			//m_spMechElemDisc->contact_forces_elem_ips_avg(locForce, elem,
			//		assoElement[0], sideCoPos, numElemCorners, locU, vActiveSetLoc);

			m_spMechElemDisc->contact_forces_elem_midpoint(locContactForce, sideElem,
					associatedElems[0], sideCoPos, locU, vActiveSetLoc);

			/*for(size_t fct=0; fct < locForce.num_all_fct(); ++fct){
				for(size_t dof=0; dof < locForce.num_all_dof(fct); ++dof)
				{
					UG_LOG("locForce(fct " << fct << ", dof " << dof <<
							"): " << locForce.value(fct,dof) << "\n");
					//const size_t globIndex = indF.index(fct, dof);
					//const size_t globComp = indF.comp(fct, dof);
					//BlockRef(force[globIndex], globComp) = locForce.value(fct,dof);
				}
			}*/

			// 	send local to global contact force
			AddLocalVector(contactForce, locContactForce);
			//	TODO: not add forces, set forces?!
		}
	} // end (elements)
}

template <typename TDomain, typename TGridFunction>
void ContactSmallStrainMechanics<TDomain, TGridFunction>::lagrange_multiplier(
		TGridFunction& contactForce, const TGridFunction& u,
		vector<DoFIndex> vActiveSetGlob,
		vector<int> vActiveSubsets)
{
	if (m_spMechElemDisc.invalid())
		UG_THROW("No element discretization set in "
					"ContactSmallStrainMechanics:lagrange_multiplier \n");

	UG_LOG("Active subsets: " << vActiveSubsets.size() << "\n");

	SmartPtr<DoFDistribution> dd = contactForce.dof_distribution();

	for (vector<int>::iterator siContact = vActiveSubsets.begin();
			siContact != vActiveSubsets.end(); ++siContact)
	{
		UG_LOG("siContact: " << *siContact << "\n");
		const int subsetDim = DimensionOfSubset(*dd->subset_handler(), *siContact);

		switch(subsetDim)
		{
		case 1:
			contact_forces_elem<RegularEdge>
				(dd->template begin<RegularEdge>(*siContact), dd->template end<RegularEdge>(*siContact),
						u, contactForce, vActiveSetGlob);
			break;
		case 2:
			contact_forces_elem<Triangle>
				(dd->template begin<Triangle>(*siContact), dd->template end<Triangle>(*siContact),
						u, contactForce, vActiveSetGlob);
			contact_forces_elem<Quadrilateral>
				(dd->template begin<Quadrilateral>(*siContact), dd->template end<Quadrilateral>(*siContact),
						u, contactForce, vActiveSetGlob);
			break;
		default:
			UG_THROW("ContactSmallStrainMechanics::lagrange_multiplier:"
				"SubsetDimension "<< subsetDim <<" (subset="<< *siContact <<") not supported.");
		}
	}

}

}//	end of namespace SmallStrainMechanics
}//	end of namespace ug


#endif /* CONTACT_IMPL_H_ */
