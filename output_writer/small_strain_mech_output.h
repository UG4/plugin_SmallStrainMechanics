/*
 * small_strain_mech_output.h
 *
 *  Created on: 21.03.2013
 *      Author: raphaelprohl
 */

/*#ifndef SMALL_STRAIN_MECH_OUTPUT_H_
#define SMALL_STRAIN_MECH_OUTPUT_H_

#include "mech_output_writer.h"

namespace ug{
namespace SmallStrainMechanics{

/// \addtogroup small_strain_mechanics
/// \{

template<typename TDomain, typename TGridFunction>
class SmallStrainMechOutput:
	public IMechOutputWriter<TDomain>
{
	private:
	///	Base class type
		typedef IMechOutputWriter<TDomain> base_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	constructor
		SmallStrainMechOutput():
			m_bIP_values_written(false), m_stressEV(false),
			m_normalStress(false){};

	///	computes the normal components of stress tensor sigma, the deviatoric part of sigma
	///	and the linearized strain tensor epsilon
		void normal_stresses_strains(TGridFunction& devSigma, TGridFunction& sigma,
				TGridFunction& epsilon, const TGridFunction& u);

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
	//	TODO: implement the 2D version by means of '=0' for the third parameter
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

	public:
		//////////////////////////
		//	INTERFACE METHODS
		//////////////////////////
		virtual void preprocess();

		template<typename TFEGeom>
		virtual void post_timestep(const number time, SmartPtr<TDomain> dom,
				TFEGeom& geo, GeometricObject* elem);

		virtual void postprocess();

	private:
		template<typename TFEGeom>
		void stress_eigenvalues_near_point(const number time, TFEGeom& geo);
		template<typename TFEGeom>
		void normal_stress_near_point(const number time, TFEGeom& geo);
		template<typename TElem>
		void normal_stress_strain_loc(LocalVector& locDevSigma, LocalVector& locSigma,
			LocalVector& locEps, TElem* elem, const LocalVector& locU);

	private:
	///	material law
		using base_type::m_spMatLaw;

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

// end group small_strain_mechanics
/// \}

} //end of namespace SmallStrainMechanics
} //end of namespace ug

#include "small_strain_mech_output_impl.h"

#endif // SMALL_STRAIN_MECH_OUTPUT_H_ */
