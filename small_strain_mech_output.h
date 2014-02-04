/*
 * small_strain_mech_output.h
 *
 *  Created on: 21.03.2013
 *      Author: raphaelprohl
 */

#ifndef SMALL_STRAIN_MECH_OUTPUT_H_
#define SMALL_STRAIN_MECH_OUTPUT_H_


namespace ug{
namespace SmallStrainMechanics{

/// \addtogroup small_strain_mechanics
/// \{

template <typename TDomain, typename TGridFunction>
class SmallStrainMechOutput
{
	public:
	///	constructor
		SmallStrainMechOutput(SmartPtr<SmallStrainMechanicsElemDisc<TDomain> > spElemDisc);

	///	computes the normal components of stress tensor sigma, the deviatoric part of sigma
	///	and the linearized strain tensor epsilon
		void normal_stresses_strains(TGridFunction& devSigma, TGridFunction& sigma,
				TGridFunction& epsilon, const TGridFunction& u);

	private:
		SmartPtr<SmallStrainMechanicsElemDisc<TDomain> > m_spElemDisc;
};

// end group small_strain_mechanics
/// \}

} //end of namespace SmallStrainMechanics
} //end of namespace ug

#include "small_strain_mech_output_impl.h"

#endif /* SMALL_STRAIN_MECH_OUTPUT_H_ */
