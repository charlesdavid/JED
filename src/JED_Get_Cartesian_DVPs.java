package jed;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import Jama.Matrix;
import bits.kde.KernelDensityEstimate2d;

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

	static String directory, out_dir, description, type, file_name_head, path, Q = "COV", R = "CORR", PC = "P_CORR";
	static int number_of_modes, reference_column, number_of_conformations;
	static int ROWS, COLS, number_of_residues;
	static Matrix FE, projections, normed_projections, weighted_normed_projections, weighted_projections, input_coords, delta_vector_series, top_evectors;
	public List<Double> eigenvalues;
	private boolean exist, success;
	KernelDensityEstimate2d KDE;

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
	 *            The type of PCA: COV or CORR (Q or R)
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
							success = (new File(out_dir)).mkdirs();
							if (!success) System.err.println("Could not create the output directory: " + out_dir);
						}
				}
			if (type.equals(R))
				{
					out_dir = directory + "JED_RESULTS_" + description + "/cPCA/CORR/";
					exist = new File(out_dir).exists();
					if (!exist)
						{
							success = (new File(out_dir)).mkdirs();
							if (!success) System.err.println("Could not create the output directory: " + out_dir);
						}
				}
			if (type.equals(PC))
				{
					out_dir = directory + "JED_RESULTS_" + description + "/cPCA/PCORR/";
					exist = new File(out_dir).exists();
					if (!exist)
						{
							success = (new File(out_dir)).mkdirs();
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
			if (number_of_modes >= 2) get_FE();
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

	/**
	 * Method to calculate the delta_G free energy from two order parameters: DVP1 and DVP2
	 */
	private void get_FE()
		{
				{
					/* Get the first 2 DVPs to use as order parameters in the deltaG free energy calculations */
					double[] order_parameter_1 = projections.getMatrix(0, number_of_conformations - 1, 0, 0).getColumnPackedCopy();
					double[] order_parameter_2 = projections.getMatrix(0, number_of_conformations - 1, 1, 1).getColumnPackedCopy();
					double[] order_parameter_1_sorted = projections.getMatrix(0, number_of_conformations - 1, 0, 0).getColumnPackedCopy();
					double[] order_parameter_2_sorted = projections.getMatrix(0, number_of_conformations - 1, 1, 1).getColumnPackedCopy();
					/* 2D KDE using Gaussian Functions */
					double[] kde_array = new double[2 * number_of_conformations];
					for (int i = 0; i < number_of_conformations; i++)
						{
							kde_array[i + i] = order_parameter_1[i];
							kde_array[i + i + 1] = order_parameter_2[i];
						}
					Arrays.sort(order_parameter_1_sorted);
					Arrays.sort(order_parameter_2_sorted);
					double op1max = order_parameter_1_sorted[number_of_conformations - 1];
					double op2max = order_parameter_2_sorted[number_of_conformations - 1];
					double op1min = order_parameter_1_sorted[0];
					double op2min = order_parameter_2_sorted[0];
					double[] bounds = { op1min, op2min, op1max, op2max };

					KDE = KernelDensityEstimate2d.compute(kde_array, 0, (number_of_conformations - 1), bounds, null, null);

					double[] probabilities = new double[number_of_conformations];
					double[] probabilities_sorted = new double[number_of_conformations];
					for (int i = 0; i < number_of_conformations; i++)
						{
							double prob = KDE.apply(order_parameter_1[i], order_parameter_2[i]);
							probabilities[i] = prob;
							probabilities_sorted[i] = prob;
						}
					Arrays.sort(probabilities_sorted);
					final double prob_max = probabilities_sorted[number_of_conformations - 1];
					final double ln_prob_max = Math.log(prob_max);
					final double KBT = (-0.600); // Units are in kcal/mol, T = 300K (room temp)
					FE = new Matrix(number_of_conformations, 3);
					for (int i = 0; i < number_of_conformations; i++)
						{
							double prob = probabilities[i];
							double ln_prob = Math.log(prob);
							double delta_G = KBT * (ln_prob - ln_prob_max);
							if (delta_G <= 0) delta_G = 0.00;// solves the -0.0000000000000 problem...
							FE.set(i, 0, order_parameter_1[i]);
							FE.set(i, 1, order_parameter_2[i]);
							FE.set(i, 2, delta_G);
						}
					path = file_name_head + "_top_2_DVPs_delta_G_" + type + ".txt";
					Matrix_IO.write_Matrix(FE, path, 12, 6);
					projections = null;
					FE = null;
				}
		}
}
