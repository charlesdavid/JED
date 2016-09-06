package jed;

import java.io.File;
import java.util.List;

import Jama.Matrix;

/**
 * JED class JED_Get_Cartesian_DVPs: Constructs the DVPs for the Cartesian subset.
 * Note: DVPs are like Principle Components, except that the reference structure is not the mean structure.
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
 * along with this program. If not, see <http://www.gnu.org/license>
 *
 * @author Dr. Charles David
 */

public class JED_Get_Cartesian_DVPs
{

	static String directory, out_dir, description, type, file_name_head, path, Q = "COV", R = "CORR", PC = "PCORR";
	static int number_of_modes, reference_column, number_of_conformations;
	static int ROWS, COLS, number_of_residues;
	static List<Double> eigenvalues;
	static Matrix projections, normed_projections, weighted_normed_projections, weighted_projections, input_coords, delta_vector_series, top_evectors;
	static long startTime, endTime, totalTime;
	boolean exist;
	boolean success;

	/**
	 * Constructor to create the delta vector series and the delta vector projections.
	 *
	 * @param data
	 *            The transformed Coordinates Matrix
	 * @param evects
	 *            The top Cartesian eigenvectors
	 * @param evals
	 *            The top Cartesian eigenvalues
	 * @param ref
	 *            The Reference Column in the Coordinates Matrix
	 * @param dir
	 *            The working directory
	 * @param des
	 *            The job description
	 * @param T
	 *            The type of PCA: COV, CORR, PCORR (Q, R, P)
	 */
	JED_Get_Cartesian_DVPs(Matrix data, Matrix evects, List<Double> evals, int ref, String dir, String des, String T)
		{

			input_coords = data;
			top_evectors = evects;
			eigenvalues = evals;
			reference_column = ref;
			directory = dir;
			description = des;
			type = T;
			if (type.equals(Q))
				{
					out_dir = directory + "JED_RESULTS_" + description + "/cPCA/COV/";
					exist = new File(out_dir).exists();
					if (!exist)
						{
							success = new File(out_dir).mkdirs();
							if (!success) System.err.println("Could not create the output directory: " + out_dir);
						}
				}
			if (type.equals(R))
				{
					out_dir = directory + "JED_RESULTS_" + description + "/cPCA/CORR/";
					exist = new File(out_dir).exists();
					if (!exist)
						{
							success = new File(out_dir).mkdirs();
							if (!success) System.err.println("Could not create the output directory: " + out_dir);
						}
				}
			if (type.equals(PC))
				{
					out_dir = directory + "JED_RESULTS_" + description + "/cPCA/PCORR/";
					exist = new File(out_dir).exists();
					if (!exist)
						{
							success = new File(out_dir).mkdirs();
							if (!success) System.err.println("Could not create the output directory: " + out_dir);
						}
				}
			number_of_modes = top_evectors.getColumnDimension();
			number_of_conformations = input_coords.getColumnDimension();
			ROWS = top_evectors.getRowDimension();
			number_of_residues = (ROWS / 3);
			delta_vector_series = new Matrix(ROWS, number_of_conformations);
			file_name_head = out_dir + "ss_" + number_of_residues;
		}

	/**
	 * Computes the Cartesian delta vector series and writes it to file.
	 */
	public void get_Cartesian_DV_Series()
		{
			Matrix ref_col = input_coords.getMatrix(0, ROWS - 1, reference_column, reference_column);
			Matrix col = new Matrix(ROWS, 1);
			for (int b = 0; b < number_of_conformations; b++)
				{
					col = input_coords.getMatrix(0, ROWS - 1, b, b);
					Matrix delta = col.minus(ref_col);
					delta_vector_series.setMatrix(0, ROWS - 1, b, b, delta);
				}
			input_coords = null;
			System.gc();
			path = directory + "JED_RESULTS_" + description + "/cPCA/ss_" + number_of_residues + "_delta_vectors.txt";
			Matrix_IO.write_Matrix(delta_vector_series, path, 9, 3);
			get_DVPs();
		}

	private void get_DVPs()
		{
			projections = new Matrix(number_of_conformations, number_of_modes);
			normed_projections = new Matrix(number_of_conformations, number_of_modes);
			weighted_projections = new Matrix(number_of_conformations, number_of_modes);
			weighted_normed_projections = new Matrix(number_of_conformations, number_of_modes);

			for (int outer = 0; outer < number_of_modes; outer++)
				{
					for (int inner = 0; inner < number_of_conformations; inner++)
						{
							int row_index_1 = 0;
							int row_index_2 = ROWS - 1;
							Matrix data1 = top_evectors.getMatrix(row_index_1, row_index_2, outer, outer);
							Matrix data2 = delta_vector_series.getMatrix(row_index_1, row_index_2, inner, inner);
							double weight = eigenvalues.get(outer);
							Matrix vector1 = Projector.get_Normed_array(data1);
							Matrix vector2 = Projector.get_Normed_array(data2);
							double dp = Projector.get_InnerProduct(data1, data2);
							double normed_dp = Projector.get_InnerProduct(vector1, vector2);
							if (dp == Double.NaN) dp = 0.000;
							if (normed_dp == Double.NaN) normed_dp = 0.000;
							double w_dp = weight * dp;
							double weighted_normed_dp = weight * normed_dp;
							projections.set(inner, outer, dp);
							normed_projections.set(inner, outer, normed_dp);
							weighted_projections.set(inner, outer, w_dp);
							weighted_normed_projections.set(inner, outer, weighted_normed_dp);
						}
				}
			path = file_name_head + "_top_" + number_of_modes + "_DVPs_" + type + ".txt";
			Matrix_IO.write_Matrix(projections, path, 9, 3);
			projections = null;
			path = file_name_head + "_top_" + number_of_modes + "_normed_DVPs_" + type + ".txt";
			Matrix_IO.write_Matrix(normed_projections, path, 9, 3);
			normed_projections = null;
			path = file_name_head + "_top_" + number_of_modes + "_weighted_DVPs_" + type + ".txt";
			Matrix_IO.write_Matrix(weighted_projections, path, 9, 3);
			weighted_projections = null;
			path = file_name_head + "_top_" + number_of_modes + "_weighted_normed_DVPs_" + type + ".txt";
			Matrix_IO.write_Matrix(weighted_normed_projections, path, 9, 3);
			weighted_normed_projections = null;

			System.gc();
		}
}
