package jed;

import java.io.File;
import java.util.List;

import Jama.Matrix;

/**
 * JED class JED_Get_Distance_Pair_DVPs: Constructs the DVPs for the Residue Distance Pairs subset.
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
 * along with this program. If not, see <http://www.gnu.org/license
 */

public class JED_Get_Distance_Pair_DVPs
{

	String directory, out_dir, description, type, file_name_head, path, Q = "COV", R = "CORR", PC = "PCORR";
	int number_of_modes, reference_column, number_of_conformations, number_of_pairs, ROWS, COLS;
	Matrix input_coords, delta_vector_series, top_evectors, projections, normed_projections, weighted_projections, weighted_normed_projections;
	Projector P;
	List<Double> top_eigenvals;
	boolean exist, success;
	KernelDensityEstimate2d KDE;

	/* **************************************** CONSTRUCTOR ************************************************************************ */

	/**
	 * Constructor to calculate the delta vectors and the DVPs for the Distance Pair analysis.
	 *
	 * @param data
	 *            The distances for the residue pairs
	 * @param evects
	 *            The eigenvectors from the Distance Pair PCA
	 * @param top_eigenvalues
	 *            The top eigenvalues form the Distance Pair PCA
	 * @param ref
	 *            The reference column to use to create the delta vectors
	 * @param dir
	 *            The working directory
	 * @param des
	 *            The job description
	 * @param T
	 *            The model type for the PCA: COV or CORR (Q or R)
	 * @param num_modes
	 *            The number of PCA modes to process
	 */
	public JED_Get_Distance_Pair_DVPs(Matrix data, Matrix evects, List<Double> top_eigenvalues, int ref, String dir, String des, String T, int num_modes)
		{

			input_coords = data;
			top_evectors = evects;
			top_eigenvals = top_eigenvalues;
			reference_column = ref;
			directory = dir;
			description = des;
			type = T;
			if (type.equals(Q))
				{
					out_dir = directory + "JED_RESULTS_" + description + "/dpPCA/COV/";
					exist = new File(out_dir).exists();
					if (!exist) success = (new File(out_dir)).mkdirs();
				}
			if (type.equals(R))
				{
					out_dir = directory + "JED_RESULTS_" + description + "/dpPCA/CORR/";
					exist = new File(out_dir).exists();
					if (!exist) success = (new File(out_dir)).mkdirs();
				}
			if (type.equals(PC))
				{
					out_dir = directory + "JED_RESULTS_" + description + "/dpPCA/PCORR/";
					exist = new File(out_dir).exists();
					if (!exist) success = (new File(out_dir)).mkdirs();
				}
			number_of_modes = num_modes;
			ROWS = input_coords.getRowDimension();
			COLS = input_coords.getColumnDimension();
			number_of_conformations = COLS;
			delta_vector_series = new Matrix(ROWS, COLS);
			number_of_pairs = ROWS;
			file_name_head = out_dir + "ss_" + number_of_pairs + "_Residue_Pairs";
		}

	/* ******************************************************************************************************************************* */

	/**
	 * Method to calculate the matrix of delta vectors using the reference frame
	 */
	public void get_Distance_Pair_DV_Series()
		{

			Matrix ref_col = input_coords.getMatrix(0, ROWS - 1, reference_column, reference_column);
			Matrix col = new Matrix(ROWS, 1);
			for (int b = 0; b < COLS; b++)
				{
					col = input_coords.getMatrix(0, ROWS - 1, b, b);
					Matrix delta = col.minus(ref_col);
					delta_vector_series.setMatrix(0, ROWS - 1, b, b, delta);
				}
			input_coords = null;
			path = directory + "JED_RESULTS_" + description + "/dpPCA/ss_" + number_of_pairs + "_Residue_Pairs_Delta_Vectors.txt";
			Matrix_IO.write_Matrix(delta_vector_series, path, 9, 3);

			get_DVPs();
		}

	/**
	 * Method to calculate the DVPs
	 */
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
							double weight = top_eigenvals.get(outer);
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
			delta_vector_series = null;
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
