/*
 * mech_output_writer.h
 *
 *  Created on: 06.02.2014
 *      Author: raphaelprohl
 */

#ifndef MECH_OUTPUT_WRITER_H_
#define MECH_OUTPUT_WRITER_H_

namespace ug{
namespace SmallStrainMechanics{

/// \addtogroup small_strain_mechanics
/// \{

class IMechOutputWriter
{
	public:
	///	constructor
		IMechOutputWriter(){};

	///	destructor
		virtual ~IMechOutputWriter(){};

	public:
		virtual void preprocess(){};
		virtual void pre_timestep(){};
		virtual void post_timestep(){};
		virtual void postprocess(){};
};

// end group small_strain_mechanics
/// \}

} //end of namespace SmallStrainMechanics
} //end of namespace ug

#endif /* MECH_OUTPUT_WRITER_H_ */
