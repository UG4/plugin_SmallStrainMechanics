/*
 * mech_output_writer.h
 *
 *  Created on: 06.02.2014
 *      Author: raphaelprohl
 */

#ifndef MECH_OUTPUT_WRITER_H_
#define MECH_OUTPUT_WRITER_H_

using namespace std;

namespace ug{
namespace SmallStrainMechanics{

/// \addtogroup small_strain_mechanics
/// \{
/// Post-Processor (Output-Writer) for Linear Elasticity and Elasto-Plasticity problems
/**
* This class delivers some functions to write/visualize some output data like e.g.
* normal components of the stress-tensor, the maximal eigenvalue of the stress tensor
* at a given point or the equivalent plastic strain tensor.
*/

template<typename TDomain>
class MechOutputWriter
{
	private:
	///	base element type of associated domain
		typedef typename domain_traits<TDomain::dim>::grid_base_object TBaseElem;

	public:
	///	World dimension
		static const int dim = TDomain::dim;

	public:
	///	constructor
		MechOutputWriter():
			m_quadOrder(0), m_bIP_values_written(false), m_stressEV(false),
			m_normalStress(false){};

	///	get stress eigenvalues at point
		void stress_eigenvalues_at(const number CoordX,
				const number CoordY, const number CoordZ)
		{
			m_stressEV = true;
			m_evalPointEV[0] = CoordX;
			m_evalPointEV[1] = CoordY;
			if (dim == 3)
				m_evalPointEV[2] = CoordZ;
		}

	///	get normal stresses at point
		void normal_stresses_at(const number CoordX,
				const number CoordY, const number CoordZ)
		{
			//	TODO:check, if point is within the domain!
			//	if not, set m_normalStress to false!
			m_normalStress = true;
			m_evalPointNormStress[0] = CoordX;
			m_evalPointNormStress[1] = CoordY;
			if (dim == 3)
				m_evalPointNormStress[2] = CoordZ;
		}

	///	destructor
		~MechOutputWriter(){};

	public:
		void preprocess();
		void pre_timestep()
		{ m_bIP_values_written = false;}

		template<typename TFEGeom>
		void post_timestep_elem(const number time, SmartPtr<TDomain> dom,
				TFEGeom& geo, TBaseElem* elem, const LocalVector& u);

		void post_timestep(const number time);
		void postprocess();

		void material_law(SmartPtr<IMaterialLaw<TDomain> > spMatLaw)
			{m_spMatLaw = spMatLaw;}
		void quad_order(const int quadOrder)
			{m_quadOrder = quadOrder;}

		void plastIP_elem(LocalVector& locPlastIP, TBaseElem* elem,
				const LocalVector& locU, SmartPtr<TDomain> dom);

		void eqPlastStrain_elem(LocalVector& locEqPlastStrain, TBaseElem* elem, const LocalVector& locU,
				SmartPtr<TDomain> dom);

		void normal_stress_strain_elem(LocalVector& locSigma, LocalVector& locEps, LocalVector& locStressFunc,
				LocalVector& locStrainFunc, TBaseElem* elem, const LocalVector& locU, SmartPtr<TDomain> dom);

	private:
		template<typename TFEGeom>
		void next_ips_to_point(vector<size_t>& vNextIP, const MathVector<dim>& point,
				const TFEGeom& geo);

		template<typename TFEGeom>
		void stress_eigenvalues_near_point(const number time,
				TFEGeom& geo, const LocalVector& u);

		template<typename TFEGeom>
		void normal_stress_near_point(const number time,
				TFEGeom& geo, const LocalVector& u);

	private:
	///	material law
		SmartPtr<IMaterialLaw<TDomain> > m_spMatLaw;

	///	quadrature rule
		int m_quadOrder;

	///	flag indicating, if ip-values are already written
		bool m_bIP_values_written;

	///	flags for stress-eigenvalues / normalStress at specific corner
		bool m_stressEV;
		bool m_normalStress;

	/// points, where to get stress-eigenvalues/normalStress
		MathVector<dim> m_evalPointEV;
		MathVector<dim> m_evalPointNormStress;

	///	output-file for stress eigenvalues
		FILE* m_fileStressEV;

};

//	some global output functions

template<typename TGridFunction>
void normal_stresses_strains(MechOutputWriter<typename TGridFunction::domain_type>& mechOut,
		TGridFunction& sigma, TGridFunction& epsilon, TGridFunction& stressFunc,
		TGridFunction& strainFunc, TGridFunction& u);

template <typename TGridFunction>
void plast_ip(TGridFunction& plastIP,
		MechOutputWriter<typename TGridFunction::domain_type>& mechOut,
		TGridFunction& u);

template <typename TGridFunction>
void equiv_plast_strain(TGridFunction& eqPlastStrain,
		MechOutputWriter<typename TGridFunction::domain_type>& mechOut,
		TGridFunction& u);

// end group small_strain_mechanics
/// \}

} //end of namespace SmallStrainMechanics
} //end of namespace ug

#include "mech_output_writer_impl.h"

#endif /* MECH_OUTPUT_WRITER_H_ */
