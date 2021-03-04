/*
 * Copyright (c) 2014-2020:  G-CSC, Goethe University Frankfurt
 * Author: Lukas Larisch
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

#ifndef VOIGTIAN_H_
#define VOIGTIAN_H_

namespace ug{
namespace SmallStrainMechanics{

template <typename TDomain>
class VoigtianMatrix{
public:
	static const int dim = IMaterialLaw<TDomain>::dim;

	VoigtianMatrix(){}
	
	void set_orthotropic(const number C11, const number C12, const number C13,
						 const number C22, const number C23, const number C33,
						 const number C44, const number C55, const number C66)
	{
		mat(0, 0) = C11;	mat(0, 1) = C12;	mat(0, 2) = C13;
		mat(1, 0) = C12;	mat(1, 1) = C22;	mat(1, 2) = C23;
		mat(2, 0) = C13;	mat(2, 1) = C23;	mat(2, 2) = C33;

		c44 = C44;
		c55 = C55;
		c66 = C66;
	}

	//counterclock-wise
	void rotate_z_90_deg(){
		DenseMatrix<FixedArray2<double, 3, 3> > rot;

		rot(0, 0) =  0; //cos 90
		rot(0, 1) = -1; //-sin 90
		rot(1, 0) =  1; //sin 90
		rot(1, 1) =  0; //cos 90
		rot(2, 2) =  1; //fixed axis

		mat = rot * mat;

		c44 = -c44;
		number tmp = c55;
		c55 = c66;
		c66 = tmp;
	}

	void copy_to_tensor(MathTensor4<dim, dim, dim, dim> &elastTensorFunct)
	{
		//  setze alle wWerte auf 0
		for (size_t i = 0; i < (size_t) dim; ++i)
			for (size_t j = 0; j < (size_t) dim; ++j)
				for (size_t k = 0; k < (size_t) dim; ++k)
					for (size_t l = 0; l < (size_t) dim; ++l)
						elastTensorFunct[i][j][k][l] = 0.0;

		elastTensorFunct[0][0][0][0] = mat(0, 0);
		elastTensorFunct[0][0][1][1] = mat(0, 1);
		elastTensorFunct[0][0][2][2] = mat(0, 2);
		//elastTensorFunct[0][0][1][2] = mat(0, 3);
		//elastTensorFunct[0][0][2][0] = mat(0, 4);
		//elastTensorFunct[0][0][0][1] = mat(0, 5);

		elastTensorFunct[1][1][0][0] = mat(1, 0);
		elastTensorFunct[1][1][1][1] = mat(1, 1);
		elastTensorFunct[1][1][2][2] = mat(1, 2);
		//elastTensorFunct[1][1][1][2] = mat(1, 3);
		//elastTensorFunct[1][1][2][0] = mat(1, 4);
		//elastTensorFunct[1][1][0][1] = mat(1, 5);

		elastTensorFunct[2][2][0][0] = mat(2, 0);
		elastTensorFunct[2][2][1][1] = mat(2, 1);
		elastTensorFunct[2][2][2][2] = mat(2, 2);
		//elastTensorFunct[2][2][1][2] = mat(2, 3);
		//elastTensorFunct[2][2][2][0] = mat(2, 4);
		//elastTensorFunct[2][2][0][1] = mat(2, 5);

		//elastTensorFunct[1][2][0][0] = mat(3, 0);
		//elastTensorFunct[1][2][1][1] = mat(3, 1);
		//elastTensorFunct[1][2][2][2] = mat(3, 2);
		////elastTensorFunct[1][2][1][2] = mat(3, 3);
		elastTensorFunct[1][2][1][2] = c44;
		//elastTensorFunct[1][2][2][0] = mat(3, 4);
		//elastTensorFunct[1][2][0][1] = mat(3, 5);

		//elastTensorFunct[2][0][0][0] = mat(4, 0);
		//elastTensorFunct[2][0][1][1] = mat(4, 1);
		//elastTensorFunct[2][0][2][2] = mat(4, 2);
		//elastTensorFunct[2][0][1][2] = mat(4, 3);
		////elastTensorFunct[2][0][2][0] = mat(4, 4);
		elastTensorFunct[2][0][2][0] = c55;
		//elastTensorFunct[2][0][0][1] = mat(4, 5);

		//elastTensorFunct[0][1][0][0] = mat(5, 0);
		//elastTensorFunct[0][1][1][1] = mat(5, 1);
		//elastTensorFunct[0][1][2][2] = mat(5, 2);
		//elastTensorFunct[0][1][1][2] = mat(5, 3);
		//elastTensorFunct[0][1][2][0] = mat(5, 4);
		////elastTensorFunct[0][1][0][1] = mat(5, 5);
		elastTensorFunct[0][1][0][1] = c66;

		//symmetries: c_ijkl = c_jikl = c_ijlk
		for(size_t i = 0; i < dim; ++i){
				for(size_t j = 0; j < dim; ++j){
						for(size_t k = 0; k < dim; ++k){
								for(size_t l = 0; l < dim; ++l){
									if(elastTensorFunct[i][j][k][l] != 0){
										elastTensorFunct[j][i][k][l] = elastTensorFunct[i][j][k][l];
										elastTensorFunct[j][i][l][k] = elastTensorFunct[i][j][k][l];
										elastTensorFunct[i][j][l][k] = elastTensorFunct[i][j][k][l];
									} 
								}
						}
				}
		}
	}

private:
	DenseMatrix<FixedArray2<double, dim, dim> > mat;
	number c44; //1/G_xy
	number c55; //1/G_yz
	number c66; //1/G_xz
};


}//	end of namespace SmallStrainMechanics
}//	end of namespace ug

#endif /* Voigtian_H_ */
