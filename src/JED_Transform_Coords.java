package jed;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 * JED class JED_Transform_Coords: Core class for calculating the rotation to align a frame to the reference frame.
 * Copyright (C) 2012 Dr. Charles David
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/license>.
 * 
 * @author Dr. Charles David
 */

public class JED_Transform_Coords
{
	final int COLS = 1;
	int ROWS, number_of_alpha_carbons;
	double ref_x_centroid, ref_y_centroid, ref_z_centroid, x_centroid, y_centroid, z_centroid;
	Matrix reference_conformation, conformation, reference_centered_coordinates, centered_coordinates, transformed_vectors, eigen_values, eigen_vectors;

	/* ************************************** CONSTRUCTORS ******************************************************************************** */

	/**
	 * Constructor for aligning frames to a reference frame
	 * 
	 * @param ref_conf
	 *            The reference conformation
	 * @param conf
	 *            The conformation to transform
	 */
	JED_Transform_Coords(Matrix ref_conf, Matrix conf)
		{

			reference_conformation = ref_conf;
			conformation = conf;
			ROWS = reference_conformation.getRowDimension();
			number_of_alpha_carbons = (ROWS / 3);
			Center_Coordinates ref = new Center_Coordinates(reference_conformation);
			reference_centered_coordinates = ref.get_centered_coordinates();
			Center_Coordinates cf = new Center_Coordinates(conformation);
			centered_coordinates = cf.get_centered_coordinates();
			ref_x_centroid = ref.x_centroid;
			ref_y_centroid = ref.y_centroid;
			ref_z_centroid = ref.z_centroid;
			x_centroid = cf.x_centroid;
			y_centroid = cf.y_centroid;
			z_centroid = cf.z_centroid;
		}

	/* ************************************** METHODS ******************************************************************************** */

	/**
	 * @return The matrix of transformed coordinates. This matrix has the same dimensions as the untransformed coordinates matrix.
	 */
	public Matrix get_transformed_coords()
		{

			Matrix N = get_matrix_N();
			EigenvalueDecomposition evd = new EigenvalueDecomposition(N);
			eigen_vectors = evd.getV();
			Matrix matrix_Q = eigen_vectors.getMatrix(0, 3, 3, 3);
			Matrix rotation_matrix = get_euler(matrix_Q);
			transformed_vectors = new Matrix(ROWS, COLS);
			for (int i = 0; i < number_of_alpha_carbons; i++)
				{
					Matrix b = new Matrix(3, 1);
					b.set(0, 0, centered_coordinates.get(i, 0));
					b.set(1, 0, centered_coordinates.get(i + number_of_alpha_carbons, 0));
					b.set(2, 0, centered_coordinates.get(i + 2 * number_of_alpha_carbons, 0));

					Matrix rotation = rotation_matrix.inverse();
					Matrix r_transformed = rotation.times(b).transpose();

					transformed_vectors.set(i, 0, r_transformed.get(0, 0));
					transformed_vectors.set(i + number_of_alpha_carbons, 0, r_transformed.get(0, 1));
					transformed_vectors.set(i + 2 * number_of_alpha_carbons, 0, r_transformed.get(0, 2));
				}
			return transformed_vectors;
		}

	private Matrix get_matrix_N()
		{

			double Sxx = 0;
			double Sxy = 0;
			double Sxz = 0;
			double Syx = 0;
			double Syy = 0;
			double Syz = 0;
			double Szx = 0;
			double Szy = 0;
			double Szz = 0;

			for (int i = 0; i < number_of_alpha_carbons; i++)
				{
					double xa = reference_centered_coordinates.get(i, 0);
					double ya = reference_centered_coordinates.get(i + number_of_alpha_carbons, 0);
					double za = reference_centered_coordinates.get(i + 2 * number_of_alpha_carbons, 0);

					double xb = centered_coordinates.get(i, 0);
					double yb = centered_coordinates.get(i + number_of_alpha_carbons, 0);
					double zb = centered_coordinates.get(i + 2 * number_of_alpha_carbons, 0);

					double xx = xa * xb;
					double xy = xa * yb;
					double xz = xa * zb;
					double yx = ya * xb;
					double yy = ya * yb;
					double yz = ya * zb;
					double zx = za * xb;
					double zy = za * yb;
					double zz = za * zb;

					Sxx += xx;
					Sxy += xy;
					Sxz += xz;
					Syx += yx;
					Syy += yy;
					Syz += yz;
					Szx += zx;
					Szy += zy;
					Szz += zz;
				}

			Matrix matrix_N = new Matrix(4, 4);

			matrix_N.set(0, 0, (Sxx + Syy + Szz));
			matrix_N.set(0, 1, (Syz - Szy));
			matrix_N.set(0, 2, (Szx - Sxz));
			matrix_N.set(0, 3, (Sxy - Syx));
			matrix_N.set(1, 0, (Syz - Szy));
			matrix_N.set(1, 1, (Sxx - Syy - Szz));
			matrix_N.set(1, 2, (Sxy + Syx));
			matrix_N.set(1, 3, (Szx + Sxz));
			matrix_N.set(2, 0, (Szx - Sxz));
			matrix_N.set(2, 1, (Sxy + Syx));
			matrix_N.set(2, 2, (-Sxx + Syy - Szz));
			matrix_N.set(2, 3, (Syz + Szy));
			matrix_N.set(3, 0, (Sxy - Syx));
			matrix_N.set(3, 1, (Szx + Sxz));
			matrix_N.set(3, 2, (Syz + Szy));
			matrix_N.set(3, 3, (-Sxx - Syy + Szz));

			return matrix_N;
		}

	private Matrix get_euler(Matrix Q)
		{

			Matrix euler = new Matrix(3, 3);

			euler.set(0, 0, (Math.pow(Q.get(0, 0), 2) + Math.pow(Q.get(1, 0), 2) - Math.pow(Q.get(2, 0), 2) - Math.pow(Q.get(3, 0), 2)));
			euler.set(0, 1, 2 * (Q.get(1, 0) * Q.get(2, 0) - Q.get(0, 0) * Q.get(3, 0)));
			euler.set(0, 2, 2 * (Q.get(1, 0) * Q.get(3, 0) + Q.get(0, 0) * Q.get(2, 0)));
			euler.set(1, 0, 2 * (Q.get(2, 0) * Q.get(1, 0) + Q.get(0, 0) * Q.get(3, 0)));
			euler.set(1, 1, (Math.pow(Q.get(0, 0), 2) - Math.pow(Q.get(1, 0), 2) + Math.pow(Q.get(2, 0), 2) - Math.pow(Q.get(3, 0), 2)));
			euler.set(1, 2, 2 * (Q.get(2, 0) * Q.get(3, 0) - Q.get(0, 0) * Q.get(1, 0)));
			euler.set(2, 0, 2 * (Q.get(3, 0) * Q.get(1, 0) - Q.get(0, 0) * Q.get(2, 0)));
			euler.set(2, 1, 2 * (Q.get(3, 0) * Q.get(2, 0) + Q.get(0, 0) * Q.get(1, 0)));
			euler.set(2, 2, (Math.pow(Q.get(0, 0), 2) - Math.pow(Q.get(1, 0), 2) - Math.pow(Q.get(2, 0), 2) + Math.pow(Q.get(3, 0), 2)));

			return euler;
		}
}
