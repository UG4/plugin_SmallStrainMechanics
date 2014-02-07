/*
 * small_strain_mech_tools_impl.h
 *
 *  Created on: 16.05.2012
 *      Author: raphaelprohl, Andreas Vogel
 */

#include "common/util/string_table_stream.h"

#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

using namespace std;

namespace ug {
namespace SmallStrainMechanics{

/// \addtogroup small_strain_mechanics
/// \{

template<typename TDomain>
void
SmallStrainMechanicsElemDisc<TDomain>::
print_mat_constants(const number lambda, const number mu,
		const number E, const number v)
{
	UG_LOG("\n Material parameter: \n");
	UG_LOG( "1. Lame constant lambda: " << lambda << "\n");
	UG_LOG( "2. Lame constant mue (sometimes 'G', shear modulus): " << mu << "\n");
	UG_LOG( "This setting equals: \n");
	UG_LOG( "young modulus E: " << E << "\n");
	UG_LOG( "poisson ratio v: " << v << "\n \n");
}

/*template<typename TDomain>
number
SmallStrainMechanicsElemDisc<TDomain>::
MatDeviatorTrace(const MathMatrix<dim, dim>& mat, MathMatrix<dim, dim>& dev)
{
	number trace = Trace(mat);

	//	compute the deviatoric part of mat
	for (size_t i = 0; i < (size_t) dim; ++i){
		for (size_t j = 0; j < (size_t) dim; ++j){
			dev[i][j] = mat[i][j];
		}
		dev[i][i] -= 1.0 / dim * trace;
	}

	return trace;
}*/

// end group small_strain_mechanics
/// \}

} //end of namespace SmallStrainMechanics
} //end of namespace ug
