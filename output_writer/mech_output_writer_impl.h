/*
 * mech_output_writer_impl.h
 *
 *  Created on: 07.02.2014
 *      Author: raphaelprohl
 */

#ifndef MECH_OUTPUT_WRITER_IMPL_H_
#define MECH_OUTPUT_WRITER_IMPL_H_

#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

// module intern headers
#include "mech_output_writer.h"

namespace ug {
namespace SmallStrainMechanics{

template <typename TGridFunction>
void plast_ip(TGridFunction& plastIP,
		MechOutputWriter<typename TGridFunction::domain_type>& mechOut,
		TGridFunction& u)
{
	static const int dim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;

	// 	local indices and local algebra
	LocalIndices ind; LocalVector locU, locPlastIP;

	const_iterator iter = u.template begin<grid_base_object>();
	const_iterator end = u.template end<grid_base_object>();

	//	loop over all elements
	for(;iter != end; ++iter)
	{
		//	get element
		grid_base_object* elem = *iter;

		// 	get global indices
		u.indices(elem, ind);

		// 	adapt local algebra
		locU.resize(ind); locPlastIP.resize(ind);

		//	reset contribution of this element
		locPlastIP = 0.0;

		//	local vector extract -> locU
		GetLocalVector(locU, u); GetLocalVector(locPlastIP, plastIP);

		mechOut.plastIP_elem(locPlastIP, elem, locU,
				u.domain());

		// 	send local to global plastIP
		AddLocalVector(plastIP, locPlastIP);
	}
}

template <typename TGridFunction>
void equiv_plast_strain(TGridFunction& eqPlastStrain,
		MechOutputWriter<typename TGridFunction::domain_type>& mechOut,
		TGridFunction& u)
{
	static const int dim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;

	// 	local indices and local algebra
	LocalIndices ind; LocalVector locU, locEqPlastStrain;

	const_iterator iter = u.template begin<grid_base_object>();
	const_iterator end = u.template end<grid_base_object>();

	//	loop over all elements
	for(;iter != end; ++iter)
	{
		//	get element
		grid_base_object* elem = *iter;

		// 	get global indices
		u.indices(elem, ind);

		// 	adapt local algebra
		locU.resize(ind); locEqPlastStrain.resize(ind);

		//	reset contribution of this element
		locEqPlastStrain = 0.0;

		//	local vector extract -> locU
		GetLocalVector(locU, u); GetLocalVector(locEqPlastStrain, eqPlastStrain);

		mechOut.eqPlastStrain_elem(locEqPlastStrain, elem, locU, u.domain());

		// 	send local to global eqPlastStrain
		AddLocalVector(eqPlastStrain, locEqPlastStrain);
	}
}

template<typename TGridFunction>
void
normal_stresses_strains(MechOutputWriter<typename TGridFunction::domain_type>& mechOut,
		TGridFunction& sigma, TGridFunction& epsilon,
		TGridFunction& stressFunc, TGridFunction& strainFunc, TGridFunction& u)
{
	static const int dim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;

	// 	local indices and local algebra
	LocalIndices indU, indEps, indSig, indStressFunc, indStrainFunc;
	LocalVector locU, locSig, locEps, locStressFunc, locStrainFunc;

	const_iterator iter = u.template begin<grid_base_object>();
	const_iterator end = u.template end<grid_base_object>();

	//	loop over all elements
	for(;iter != end; ++iter)
	{
		//	get element
		grid_base_object* elem = *iter;

		// 	get global indices
		u.indices(elem, indU); epsilon.indices(elem, indEps);
		sigma.indices(elem, indSig); stressFunc.indices(elem, indStressFunc);
		strainFunc.indices(elem, indStrainFunc);

		// 	adapt local algebra
		locU.resize(indU); locEps.resize(indEps);
		locSig.resize(indSig); locStressFunc.resize(indStressFunc);
		locStrainFunc.resize(indStrainFunc);

		//	reset contribution of this element
		locStressFunc = 0.0; locStrainFunc = 0.0; locSig = 0.0; locEps = 0.0;

		//	local vector extract -> locU
		GetLocalVector(locU, u);

		//	call method of solid mechanics output writer
		mechOut.normal_stress_strain_elem(locSig, locEps, locStressFunc, locStrainFunc, elem,
				locU, u.domain());

		// 	send local to global
		//	sigma (cauchy-stress tensor), epsilon (linearized strain tensor), devSigma (deviatoric part of Sigma)
		AddLocalVector(stressFunc, locStressFunc);
		AddLocalVector(strainFunc, locStrainFunc);
		AddLocalVector(sigma, locSig);
		AddLocalVector(epsilon, locEps);
	}
}

template <typename TGridFunction>
void invariants_kirchhoff_stress(TGridFunction& invarKirchhoffStress,
		MechOutputWriter<typename TGridFunction::domain_type>& mechOut,
		TGridFunction& u)
{
	static const int dim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;

	// 	local indices and local algebra
	LocalIndices ind, indInvarKirchhoffStress;
	LocalVector locU, locInvarKirchhoffStress;

	const_iterator iter = u.template begin<grid_base_object>();
	const_iterator end = u.template end<grid_base_object>();

	//	loop over all elements
	for(;iter != end; ++iter)
	{
		//	get element
		grid_base_object* elem = *iter;

		// 	get global indices
		u.indices(elem, ind); invarKirchhoffStress.indices(elem, indInvarKirchhoffStress);

		// 	adapt local algebra
		locU.resize(ind); locInvarKirchhoffStress.resize(indInvarKirchhoffStress);

		//	local vector extract -> locU
		GetLocalVector(locU, u); GetLocalVector(locInvarKirchhoffStress, invarKirchhoffStress);

		//	reset contribution of this element
		locInvarKirchhoffStress = 0.0;

		mechOut.invariants_kirchhoff_stress_elem(locInvarKirchhoffStress, elem, locU,
				u.domain());

		// 	send local to global tensor
		AddLocalVector(invarKirchhoffStress, locInvarKirchhoffStress);
	}
}

///	count how many integration points are plastified
template<typename TDomain>
void
MechOutputWriter<TDomain>::
plastIP_elem(LocalVector& locPlastIP, TBaseElem* elem,
		const LocalVector& locU, SmartPtr<TDomain> dom)
{
	//	get vertices and extract corner coordinates
	typedef typename TDomain::position_accessor_type position_accessor_type;
	position_accessor_type& aaPos = dom->position_accessor();
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	const size_t numVertices = elem->num_vertices();
	for (size_t i = 0; i < numVertices; ++i) {
		coCoord[i] = aaPos[elem->vertex(i)];
	};

	//	update the geometry for the Element elem
	DimFEGeometry<dim> geo;
   	try{
		geo.update(elem, &(coCoord[0]),
				LFEID(LFEID::LAGRANGE, dim, 1), m_quadOrder);
	}
   	UG_CATCH_THROW("SmallStrainMechOutput::plastIP_elem:"
					" Cannot update Finite Element Geometry.");

	//	get all neighbor elems which share a vertex with the given element 'elem'
	typename TDomain::grid_type& grid = *(dom->grid());
	typedef typename vector<TBaseElem*>::iterator neighborElemIter;
	vector<TBaseElem*> vNeighborElems;
	CollectNeighbors(vNeighborElems, elem, grid, NHT_VERTEX_NEIGHBORS);

	//  pointer to internal variable of current elem
	m_spMatLaw->internal_vars(elem);

	/*MathMatrix<dim, dim> GradU;
	number plasticIP = 0.0;
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//	get displacementGradient (GradU)
		m_spMatLaw->template DisplacementGradient<DimFEGeometry<dim> >(GradU, ip, geo, locU);

		number gamma = m_spMatLaw->plastic_multiplier(ip, GradU);
		if (gamma > 0.0)
			plasticIP += 1.0;
		if (gamma < 0.0)
			UG_THROW("gamma: " << gamma << "in plastIP_elem \n");
	}*/

	MathMatrix<dim, dim> plastStrain;
	number plasticIP = 0.0;
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		plastStrain = *(m_spMatLaw->inelastic_strain_tensor(ip));

		number plastStrainNorm = MatFrobeniusNorm(plastStrain);
		if (plastStrainNorm > 0.0)
			plasticIP += 1.0;
	}

	for (size_t co = 0; co < elem->num_vertices(); ++co)
	{
		//	init 'elemsWithCo' with 1, because co lies in 'elem'
		size_t elemsWithCo = 1;

		//	iterate neighbor elems and count how many elems contain the corner 'co'
		for (neighborElemIter it = vNeighborElems.begin();
				it != vNeighborElems.end(); ++it){
			if (ContainsPoint(*it, coCoord[co], aaPos) == true)
				elemsWithCo++;
		}

		//	scaling factor for averaging values out of all ips
		//	of the associated elements of corner co
		size_t scaleFac = elemsWithCo * geo.num_ip();

		locPlastIP(0,co) = plasticIP / scaleFac;
	}
}

///	compute the equivalent plastic strain elementwise (see Simo/Hughes p. 90 for the definition)
template<typename TDomain>
void
MechOutputWriter<TDomain>::
eqPlastStrain_elem(LocalVector& locEqPlastStrain, TBaseElem* elem,
		const LocalVector& locU, SmartPtr<TDomain> dom)
{
	//	get vertices and extract corner coordinates
	typedef typename TDomain::position_accessor_type position_accessor_type;
	position_accessor_type& aaPos = dom->position_accessor();
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	const size_t numVertices = elem->num_vertices();
	for (size_t i = 0; i < numVertices; ++i) {
		coCoord[i] = aaPos[elem->vertex(i)];
	};

	//	prepare geometry for type and order
	DimFEGeometry<dim> geo;
   	try{
		geo.update(elem, &(coCoord[0]),
				LFEID(LFEID::LAGRANGE, dim, 1), m_quadOrder);
	}
   	UG_CATCH_THROW("SmallStrainMechOutput::eqPlastStrain_elem:"
					" Cannot update Finite Element Geometry.");

	//	get all neighbor elems which share a vertex with the given element 'elem'
	typename TDomain::grid_type& grid = *(dom->grid());
	typedef typename vector<TBaseElem*>::iterator neighborElemIter;
	vector<TBaseElem*> vNeighborElems;
	CollectNeighbors(vNeighborElems, elem, grid, NHT_VERTEX_NEIGHBORS);

	//  pointer to internal variable of current elem
	m_spMatLaw->internal_vars(elem);

	MathMatrix<dim, dim> GradU;
	number SumEqPlastStrainIP = 0.0;
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//	get displacementGradient (GradU)
		m_spMatLaw->template DisplacementGradient<DimFEGeometry<dim> >(GradU, ip, geo, locU);

		number gamma = m_spMatLaw->plastic_multiplier(ip, GradU);

		//	sum over all ip`s
		SumEqPlastStrainIP += gamma * sqrt(2.0/3.0);
	}

	for (size_t co = 0; co < elem->num_vertices(); ++co)
	{
		//	init 'elemsWithCo' with 1, because co lies in 'elem'
		size_t elemsWithCo = 1;

		//	iterate neighbor elems and count how many elems contain the corner 'co'
		for (neighborElemIter it = vNeighborElems.begin();
				it != vNeighborElems.end(); ++it){
			if (ContainsPoint(*it, coCoord[co], aaPos) == true)
				elemsWithCo++;
		}

		//	scaling factor for averaging values out of all ips
		//	of the associated elements of corner co
		size_t scaleFac = elemsWithCo * geo.num_ip();

		locEqPlastStrain(0, co) = SumEqPlastStrainIP / scaleFac;
	}
}
template<typename TDomain>
void
MechOutputWriter<TDomain>::
normal_stress_strain_elem(LocalVector& locSigma, LocalVector& locEps, LocalVector& locStressFunc,
		LocalVector& locStrainFunc, TBaseElem* elem, const LocalVector& locU,
		SmartPtr<TDomain> dom)
{
//	get all neighbor elems which share a vertex with the given element 'elem'
	typename TDomain::grid_type& grid = *(dom->grid());
	vector<TBaseElem*> vNeighborElems;
	CollectNeighbors(vNeighborElems, elem, grid, NHT_VERTEX_NEIGHBORS);
	typedef typename vector<TBaseElem*>::iterator neighborElemIter;

	typedef typename TDomain::position_accessor_type position_accessor_type;
	position_accessor_type& aaPos = dom->position_accessor();

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get vertices and extract corner coordinates
	const size_t numVertices = elem->num_vertices();
	for (size_t i = 0; i < numVertices; ++i) {
		coCoord[i] = aaPos[elem->vertex(i)];
	};

	//	prepare geometry for type and order
	DimFEGeometry<dim> geo;
   	try{
		geo.update(elem, &(coCoord[0]),
				LFEID(LFEID::LAGRANGE, dim, 1), m_quadOrder);
	}
   	UG_CATCH_THROW("SmallStrainMechOutput::normal_stress_strain_elem:"
					" Cannot update Finite Element Geometry.");

	//  pointer to internal variable of current elem
	m_spMatLaw->internal_vars(elem);

	MathMatrix<dim, dim> GradU, eps, sigma, devSigma, SumEpsIP, SumSigmaIP, inelasticStrainTens, eps_elast;
	number normDevSig, normSig, normEps, SumNormDevSigIP, SumNormSigIP, SumNormEpsIP,
	SumVolSigmaIP, normElastStrain, SumNormElastStrainIP, normInelastStrain, SumNormInelastStrainIP;

	//	init summation-values
	SumEpsIP = 0.0; SumSigmaIP = 0.0; SumNormDevSigIP = 0.0; SumNormSigIP = 0.0;
	SumNormEpsIP = 0.0; SumVolSigmaIP = 0.0; SumNormElastStrainIP = 0.0; SumNormInelastStrainIP = 0.0;

	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//	get displacementGradient (GradU)
		m_spMatLaw->template DisplacementGradient<DimFEGeometry<dim> >(GradU, ip, geo, locU);
		//	get Cauchy stress-tensor (sigma) at ip
		m_spMatLaw->stressTensor(sigma, ip, GradU);

		// 	get inelastic strain tensor at ip
		inelasticStrainTens = *(m_spMatLaw->inelastic_strain_tensor(ip));

		// 	get linearized strain tensor (eps) at ip
		for (size_t i = 0; i < (size_t) dim; ++i)
			for (size_t J = 0; J < (size_t) dim; ++J)
			{
				eps[i][J] = 0.5 * (GradU[i][J] + GradU[J][i]);
				eps_elast[i][J] = eps[i][J] - inelasticStrainTens[i][J];
			}

		//	add values to local vector
		for (size_t i = 0; i < (size_t) dim; ++i)
		{
			SumEpsIP[i][i] += eps[i][i];
			SumSigmaIP[i][i] += sigma[i][i];
		}

		//	stress functions
		number trSigma = MatDeviatorTrace(sigma, devSigma);
		normDevSig = MatFrobeniusNorm(devSigma);
		normSig = MatFrobeniusNorm(sigma);
		SumVolSigmaIP += 1.0/3.0 * trSigma;
		SumNormDevSigIP += normDevSig;
		SumNormSigIP += normSig;

		//	strain functions
		normEps = MatFrobeniusNorm(eps);
		SumNormEpsIP += normEps;

		normElastStrain = MatFrobeniusNorm(eps_elast);
		SumNormElastStrainIP += normElastStrain;

		normInelastStrain = MatFrobeniusNorm(inelasticStrainTens);
		SumNormInelastStrainIP += normInelastStrain;

	} //end (ip)

	// loop corner
	for (size_t co = 0; co < numVertices; ++co)
	{
		//	TODO: is this geometry update necessary here?

		//	update geometry
		try{
			geo.update(elem, &(coCoord[0]),
					LFEID(LFEID::LAGRANGE, dim, 1), m_quadOrder);
		}
		UG_CATCH_THROW("SmallStrainMechOutput::normal_stress_strain_loc:"
						" Cannot update Finite Element Geometry.");

		//	init 'elemsWithCo' with 1, because co lies in 'elem'
		size_t elemsWithCo = 1;

		//	iterate neighbor elems and count how many elems contain the corner 'co'
		for (neighborElemIter it = vNeighborElems.begin();
				it != vNeighborElems.end(); ++it){
			if (ContainsPoint(*it, coCoord[co], aaPos) == true)
				elemsWithCo++;
		}

		//	scaling factor for averaging values out of all ips
		//	of the associated elements of corner co
		size_t scaleFac = elemsWithCo * geo.num_ip();

		//	add values to local vector
		for (size_t i = 0; i < (size_t) dim; ++i)
		{
			locEps(i, co) += SumEpsIP[i][i] / scaleFac;
			locSigma(i, co) += SumSigmaIP[i][i] / scaleFac;
		}

		locStressFunc(0, co) += SumNormSigIP / scaleFac;
		locStressFunc(1, co) += SumNormDevSigIP / scaleFac;
		locStressFunc(2, co) += SumVolSigmaIP / scaleFac;

		locStrainFunc(0, co) += SumNormEpsIP / scaleFac;
		locStrainFunc(1, co) += SumNormElastStrainIP / scaleFac;
		locStrainFunc(2, co) += SumNormInelastStrainIP / scaleFac;
	} //end(co)


	//	temporarily output for linear-elastic-benchmark-output:
	for (size_t co = 0; co < numVertices; ++co)
	{
		if ((coCoord[co][0] == 90.0) && (coCoord[co][1] == 0.0))
		{
			UG_LOG("At 2 (90.0,0.0): \n");
			for (size_t i = 0; i < (size_t) dim; ++i)
				UG_LOG("u("<< i << "," << co << "): " << locU(i,co) << "\n");

			UG_LOG("\n");

			vector<size_t> vNextIP;
			next_ips_to_point(vNextIP, coCoord[co], geo);

			MathMatrix<dim, dim> sigma;
			for (vector<size_t>::iterator it = vNextIP.begin();
							it != vNextIP.end(); ++it)
			{
				//	get displacementGradient (GradU) at ip
				m_spMatLaw->template DisplacementGradient<DimFEGeometry<dim> >(GradU, *it, geo, locU);
				//	get the Cauchy Stress Tensor (sigma) at ip
				m_spMatLaw->stressTensor(sigma, *it, GradU);

				UG_LOG("At " << geo.global_ip(*it) << ": \n");
				UG_LOG("sigma_yy: " << sigma[1][1] << "\n");
				UG_LOG("sigma: " << sigma << "\n");
			}
		}

		if ((coCoord[co][0] == 100.0) && (coCoord[co][1] == 100.0))
		{
			UG_LOG("At 4 (100.0,100.0): \n");
			for (size_t i = 0; i < (size_t) dim; ++i)
				UG_LOG("u("<< i << "," << co << "): " << locU(i,co) << "\n");

			UG_LOG("\n");
		}

		if ((coCoord[co][0] == 0.0) && (coCoord[co][1] == 100.0))
		{
			UG_LOG("At 5 (0.0,100.0): \n");
			for (size_t i = 0; i < (size_t) dim; ++i)
				UG_LOG("u("<< i << "," << co << "): " << locU(i,co) << "\n");

			UG_LOG("\n");
		}

	};
}

///	compute the invariants of the kirchhoff stress tensor per element
template<typename TDomain>
void
MechOutputWriter<TDomain>::
invariants_kirchhoff_stress_elem(LocalVector& locInvarKirchhoffStress,
		TBaseElem* elem, const LocalVector& u, SmartPtr<TDomain> dom)
{
	//	get vertices and extract corner coordinates
	typedef typename TDomain::position_accessor_type position_accessor_type;
	position_accessor_type& aaPos = dom->position_accessor();
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	const size_t numVertices = elem->num_vertices();
	for (size_t i = 0; i < numVertices; ++i) {
		coCoord[i] = aaPos[elem->vertex(i)];
	};

	//	prepare geometry for type and order
	DimFEGeometry<dim> geo;
	try{
		geo.update(elem, &(coCoord[0]),
				LFEID(LFEID::LAGRANGE, dim, 1), m_quadOrder);
	}
	UG_CATCH_THROW("MechOutputWriterFinite::invariants_kirchhoff_stress_elem:"
					" Cannot update Finite Element Geometry.");

	//	get all neighbor elems which share a vertex with the given element 'elem'
	typename TDomain::grid_type& grid = *(dom->grid());
	typedef typename vector<TBaseElem*>::iterator neighborElemIter;
	vector<TBaseElem*> vNeighborElems;
	CollectNeighbors(vNeighborElems, elem, grid, NHT_VERTEX_NEIGHBORS);

	//  pointer to internal variable of current elem
	m_spMatLaw->internal_vars(elem);

	MathMatrix<dim, dim> tauTens, tauTens2, GradU, sigma, F, Ident;
	MatIdentity(Ident);
	for (size_t co = 0; co < elem->num_vertices(); ++co)
	{
		number invar1_ip = 0.0;
		number invar2_ip = 0.0;
		number invar3_ip = 0.0;

		//	loop ip`s and determine the next ip`s to corner co!
		vector<size_t> vNextIP;
		next_ips_to_point(vNextIP, coCoord[co], geo);

		for (vector<size_t>::iterator it = vNextIP.begin();
								it != vNextIP.end(); ++it)
		{
			//	get deformation gradient at ip
			m_spMatLaw->template DisplacementGradient<DimFEGeometry<dim> >(GradU, *it, geo, u);

			//	get Cauchy stress-tensor (sigma) at ip
			m_spMatLaw->stressTensor(sigma, *it, GradU);

			//	compute the kirchhoff stress-tensor at ip
			MatAdd(F, GradU, Ident);
			number detF = Determinant(F);
			MatScale(tauTens, 1.0/detF, sigma);


			MatMultiply(tauTens2, tauTens, tauTens);
			number traceTauTens_ip = Trace(tauTens);
			number traceTauTens2_ip = Trace(tauTens2);

			invar1_ip += traceTauTens_ip;
			invar2_ip += 0.5 * (traceTauTens_ip * traceTauTens_ip - traceTauTens2_ip);
			invar3_ip += Determinant(tauTens);
		}

		//	init 'elemsWithCo' with 1, because co lies in 'elem'
		size_t elemsWithCo = 1;

		//	iterate neighbor elems and count how many elems contain the corner 'co'
		for (neighborElemIter it = vNeighborElems.begin();
				it != vNeighborElems.end(); ++it){
			if (ContainsPoint(*it, coCoord[co], aaPos) == true)
				elemsWithCo++;
		}

		//	scaling factor for averaging values out of all ips
		//	of the associated elements of corner co
		const size_t scaleFac = elemsWithCo * vNextIP.size();

		locInvarKirchhoffStress(0, co) = invar1_ip/scaleFac;
		locInvarKirchhoffStress(1, co) = invar2_ip/scaleFac;
		locInvarKirchhoffStress(2, co) = invar3_ip/scaleFac;
	}
}

template<typename TDomain>
void
MechOutputWriter<TDomain>::
preprocess()
{
	/*if (m_stressEV)
	{
		#ifdef UG_PARALLEL
			if (pcl::ProcRank() == 0){
				m_fileStressEV = fopen("sig_eigen.dat", "w");
			}
		#else
			m_fileStressEV = fopen("sig_eigen.dat", "w");
		#endif
	}*/
}

template<typename TDomain>
template<typename TFEGeom>
void
MechOutputWriter<TDomain>::
post_timestep_elem(const number time, SmartPtr<TDomain> dom, TFEGeom& geo,
		TBaseElem* elem, const LocalVector& u)
{
	typedef typename TDomain::position_accessor_type position_accessor_type;
	position_accessor_type& aaPos = dom->position_accessor();

	if (m_stressEV && (!m_bIP_values_written))
		if (ContainsPoint(elem, m_evalPointEV, aaPos))
			stress_eigenvalues_near_point(time, geo, u);

	//	COMPUTE NORMAL STRESSES OF THE STRESSTENSOR SIGMA AT IP
	if (m_normalStress)
		if (ContainsPoint(elem, m_evalPointNormStress, aaPos))
			normal_stress_near_point(time, geo, u);
}

template<typename TDomain>
void
MechOutputWriter<TDomain>::
post_timestep(const number time)
{
	//	write some material law data to console
	m_spMatLaw->write_data_to_console(time);
}

template<typename TDomain>
void
MechOutputWriter<TDomain>::
postprocess()
{
	/*if (m_stressEV)
	{
		#ifdef UG_PARALLEL
			if (pcl::ProcRank() == 0){
				fclose(m_fileStressEV);
			}
		#else
			fclose(m_fileStressEV);
		#endif
	}*/
}

template<typename TDomain>
template<typename TFEGeom>
void
MechOutputWriter<TDomain>::
stress_eigenvalues_near_point(const number time, TFEGeom& geo,
		const LocalVector& u)
{
	vector<size_t> vNextIP;
	next_ips_to_point(vNextIP, m_evalPointEV, geo);

	MathMatrix<dim, dim> SYMMSig, Sig, GradU;
	for (vector<size_t>::iterator it = vNextIP.begin();
					it != vNextIP.end(); ++it)
	{
		m_spMatLaw->template DisplacementGradient<TFEGeom>(GradU, *it, geo, u);
		m_spMatLaw->stressTensor(Sig, *it, GradU);

		if (dim == 3)
		{
			MatMultiplyMTM(SYMMSig, Sig);

			MathMatrix<3, 3, number> SigTens;
			for (size_t a = 0; a < 3; ++a)
				for (size_t b = 0; b < 3; ++b){
					SigTens[a][b] = SYMMSig[a][b];
				}

			MathVector<3, number> evMin, evMed, evMax;
			number SiglambdaMin, SiglambdaMed, SiglambdaMax;

			CalculateEigenvalues(SigTens, SiglambdaMin,
					SiglambdaMed, SiglambdaMax, evMin, evMed,
					evMax);

			number sqrSigLambdaMin = sqrt(SiglambdaMin);
			number sqrSigLambdaMed = sqrt(SiglambdaMed);
			number sqrSigLambdaMax = sqrt(SiglambdaMax);

			UG_LOG("minimal EigenValueSigma: " << sqrSigLambdaMin << "\n");
			UG_LOG("medium EigenValueSigma: " << sqrSigLambdaMed << "\n");
			UG_LOG("maximal EigenValueSigma: " << sqrSigLambdaMax << "\n");

		/*	number absSiglambdamax = abs(sqrSigLambdaMax);

		#ifdef UG_PARALLEL
			if (pcl::ProcRank() == 0){
				fprintf(m_fileStressEV, "%lg %lg \n ", time, absSiglambdamax);
			}
		#else
			fprintf(m_fileStressEV, "%lg %lg \n ", time, absSiglambdamax);
		#endif*/

		}

		//	to be sure that data is only written once:
		m_bIP_values_written = true;
	} // end vNextIP-iteration
}

template<typename TDomain>
template<typename TFEGeom>
void
MechOutputWriter<TDomain>::
normal_stress_near_point(const number time, TFEGeom& geo, const LocalVector& u)
{
	vector<size_t> vNextIP;
	next_ips_to_point(vNextIP, m_evalPointNormStress, geo);

	MathMatrix<dim, dim> Sig, GradU;
	for (vector<size_t>::iterator it = vNextIP.begin();
					it != vNextIP.end(); ++it)
	{
		m_spMatLaw->template DisplacementGradient<TFEGeom>(GradU, *it, geo, u);
		m_spMatLaw->stressTensor(Sig, *it, GradU);

		UG_LOG("At " << geo.global_ip(*it) << ": \n");
		UG_LOG("sigma_xx: " << Sig[0][0] << "\n");
		UG_LOG("sigma_yy: " << Sig[1][1] << "\n");
		UG_LOG("sigma: " << Sig << "\n");
		UG_LOG("######### \n");

		//	to be sure that data is only written once:
		//m_bIP_values_written = true;
	} // end vNextIP-iteration
}

template<typename TDomain>
template<typename TFEGeom>
void
MechOutputWriter<TDomain>::
next_ips_to_point(vector<size_t>& vNextIP, const MathVector<dim>& point,
		const TFEGeom& geo)
{
	number dist = std::numeric_limits<number>::max();

	//	determine shortest distance from an ip to a point
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		number dist_ip = VecDistance(geo.global_ip(ip), point);
		if (dist_ip < dist)
			dist = dist_ip;
	}

	//	determine all ip`s with shortest distance to a point
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		number dist_ip = VecDistance(geo.global_ip(ip), point);
		if (dist_ip == dist)
			vNextIP.push_back(ip);
	}
}


} // namespace SmallStrainMechanics
} // namespace ug

#endif /* MECH_OUTPUT_WRITER_IMPL_H_ */
