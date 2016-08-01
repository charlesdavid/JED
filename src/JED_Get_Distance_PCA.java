package jed;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 * JED class JED_Get_Distance_PCA: Gets the COV and CORR PCA for the All-to-All Distance subset.
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

public class JED_Get_Distance_PCA
{

	String directory, out_dir_dPCA, out_dir_COV, out_dir_CORR, description, file_name_head, path;
	int ROWS_D, COLS, number_of_modes, number_of_residues;
	double trace_COV, trace_CORR, cond_COV, cond_CORR;
	List<Double> eigenvalues_COV, top_eigenvalues_COV, eigenvalues_CORR, top_eigenvalues_CORR;
	double[] dist_pca_mode_COV_min, dist_pca_mode_COV_max, dist_pca_mode_CORR_min, dist_pca_mode_CORR_max;
	Matrix distances, centered_data, COV_dist, CORR_dist, top_evectors_dist_COV, top_evectors_dist_CORR, residue_means, residue_std_devs;
	NumberFormat nf;
	RoundingMode rm;
	EigenvalueDecomposition evd;
	PCA pca;
	boolean success, exist;

	/* ****************************** CONSTRUCTOR ********************************************************************************* */

	/**
	 * Constructor to perform the All to All Distance PCA
	 * 
	 * @param dist
	 *            The matrix of distances for the subset
	 * @param dir
	 *            The working directory
	 * @param des
	 *            The job description
	 * @param modes
	 *            The number of distance PCA modes
	 * @param res
	 */
	JED_Get_Distance_PCA(Matrix dist, String dir, String des, int modes, int res)
		{

			directory = dir;
			description = des;

			out_dir_dPCA = directory + "JED_RESULTS_" + description + "/dPCA/";
			exist = new File(out_dir_dPCA).exists();
			if (!exist) success = (new File(out_dir_dPCA)).mkdirs();

			out_dir_COV = out_dir_dPCA + "/COV/";
			exist = new File(out_dir_COV).exists();
			if (!exist) success = (new File(out_dir_COV)).mkdirs();

			out_dir_CORR = out_dir_dPCA + "CORR/";
			exist = new File(out_dir_CORR).exists();
			if (!exist) success = (new File(out_dir_CORR)).mkdirs();

			distances = dist;
			number_of_residues = res;
			number_of_modes = modes;
			ROWS_D = (number_of_residues * (number_of_residues - 1) / 2);
			COLS = number_of_modes;

			nf = NumberFormat.getInstance();
			rm = RoundingMode.HALF_UP;
			nf.setRoundingMode(rm);
			nf.setMaximumFractionDigits(3);
			nf.setMinimumFractionDigits(3);
		}

	/* ******************************** DRIVER METHODS ************************************************************************************** */

	/**
	 * Gets the All to ALl Distance PCA using COV and CORR models
	 */
	public void get_Distance_PCA()
		{

			pca = new PCA(distances);

			distances = null;
			System.gc();

			Do_Cov();
			Do_Corr();
		}

	private void Do_Cov()
		{

			COV_dist = pca.get_covariance_matrix_elegant();

			residue_means = pca.getData_means();

			file_name_head = out_dir_dPCA + "ss_" + number_of_residues;
			path = file_name_head + "_Centroids_of_Variables.txt";
			Matrix_IO.write_Matrix(residue_means, path, 12, 6);

			residue_std_devs = pca.getData_sigmas();

			path = file_name_head + "_Std_Devs_of_Variables.txt";
			Matrix_IO.write_Matrix(residue_std_devs, path, 12, 6);

			file_name_head = out_dir_COV + "ss_" + number_of_residues;
			Matrix_IO.write_Matrix(COV_dist, file_name_head + "_COV_matrix.txt", 12, 6);
			evd = pca.get_eigenvalue_decomposition(COV_dist);

			get_eigenvalues_COV();
			write_top_evals_COV();
			get_top_evects_and_reverse_COV();
			construct_PCA_Modes_COV();
		}

	private void Do_Corr()
		{

			CORR_dist = pca.get_R_from_Q(COV_dist);
			COV_dist = null;

			file_name_head = out_dir_CORR + "ss_" + number_of_residues;
			Matrix_IO.write_Matrix(CORR_dist, file_name_head + "_CORR_matrix.txt", 12, 6);
			evd = pca.get_eigenvalue_decomposition(CORR_dist);
			get_eigenvalues_CORR();
			write_top_evals_CORR();
			get_top_evects_and_reverse_CORR();
			construct_PCA_Modes_CORR();
		}

	/* ********************************* COV METHODS ********************************************************************************** */

	private void get_eigenvalues_COV()
		{

			double[] ss_evals = evd.getRealEigenvalues();
			trace_COV = 0.000;
			eigenvalues_COV = new ArrayList<Double>();
			for (double k : ss_evals)
				{
					eigenvalues_COV.add(k);
					trace_COV += k;
				}
			Collections.sort(eigenvalues_COV, Collections.reverseOrder());
			double max = eigenvalues_COV.get(0);
			double min = eigenvalues_COV.get(ROWS_D - 1);
			cond_COV = Math.abs(max / min);
			file_name_head = out_dir_COV + "ss_" + number_of_residues;
			List_IO.write_Double_List(eigenvalues_COV, file_name_head + "_distance_eigenvalues_COV.txt", 6);

		}

	private void write_top_evals_COV()
		{
			try
				{
					file_name_head = out_dir_COV + "ss_" + number_of_residues;
					File top_ss_evals_cov = new File(file_name_head + "_top_" + number_of_modes + "_distance_eigenvalues_COV.txt");
					BufferedWriter top_ss_evals_writer = new BufferedWriter(new FileWriter(top_ss_evals_cov));
					top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", "Eigenvalue", "% Variance", "Cumulative Variance"));
					top_eigenvalues_COV = new ArrayList<Double>();
					double cumulative_variance = 0;
					for (int i = 0; i < number_of_modes; i++)
						{
							double val = eigenvalues_COV.get(i);
							double normed_val = (val / trace_COV);
							cumulative_variance += normed_val;
							top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", nf.format(val), nf.format(normed_val), nf.format(cumulative_variance)));
							top_eigenvalues_COV.add(val);
						}
					top_ss_evals_writer.close();

				} catch (IOException io)
				{
					System.err.println("IOException thrown. Could not write the file: " + file_name_head + "_top_" + number_of_modes + "_distance_eigenvalues_COV.txt");
					io.getMessage();
					io.getStackTrace();
				}
		}

	private void get_top_evects_and_reverse_COV()
		{

			Matrix ss_evectors = evd.getV();
			Matrix D = evd.getD();
			Matrix precision = ss_evectors.times(D.inverse()).times(ss_evectors.transpose());
			evd = null;
			System.gc();

			top_evectors_dist_COV = ss_evectors.getMatrix(0, ROWS_D - 1, ROWS_D - number_of_modes, ROWS_D - 1);
			Matrix modes_reversed = new Matrix(ROWS_D, COLS);
			for (int r = 0; r < COLS; r++)
				{
					Matrix col = top_evectors_dist_COV.getMatrix(0, ROWS_D - 1, COLS - 1 - r, COLS - 1 - r);
					modes_reversed.setMatrix(0, ROWS_D - 1, r, r, col);
				}
			top_evectors_dist_COV = modes_reversed;

			file_name_head = out_dir_COV + "ss_" + number_of_residues;
			path = file_name_head + "_top_" + number_of_modes + "_distance_eigenvectors_COV.txt";
			Matrix_IO.write_Matrix(top_evectors_dist_COV, path, 12, 6);
			path = file_name_head + "_PRECISION_matrix.txt";
			Matrix_IO.write_Matrix(precision, path, 12, 3);
			ss_evectors = null;
			System.gc();
		}

	private void construct_PCA_Modes_COV()
		{

			Matrix SS_pca_modes = new Matrix(ROWS_D, number_of_modes);
			Matrix SS_pca_square_modes = new Matrix(ROWS_D, number_of_modes);
			Matrix SS_weighted_pca_modes = new Matrix(ROWS_D, number_of_modes);
			Matrix SS_weighted_pca_square_modes = new Matrix(ROWS_D, number_of_modes);
			dist_pca_mode_COV_max = new double[number_of_modes];
			dist_pca_mode_COV_min = new double[number_of_modes];

			for (int a = 0; a < number_of_modes; a++)
				{
					double max = 0;
					double min = 1;

					for (int b = 0; b < ROWS_D; b++)
						{
							double d = top_evectors_dist_COV.get(b, a);
							double d_abs = Math.abs(d);
							double sq = (d * d);
							double value = top_eigenvalues_COV.get(a);
							if (value < 0) value = 0;
							double sqrt_val = Math.sqrt(value);
							double wm = sqrt_val * d_abs;
							double w_sq = value * sq;
							SS_pca_modes.set(b, a, d_abs);
							SS_weighted_pca_modes.set(b, a, wm);
							SS_pca_square_modes.set(b, a, sq);
							SS_weighted_pca_square_modes.set(b, a, w_sq);
							if (d_abs >= max)
								{
									max = d_abs;
								}
							if (d_abs <= min)
								{
									min = d_abs;
								}
						}

					dist_pca_mode_COV_max[a] = max;
					dist_pca_mode_COV_min[a] = min;
				}
			file_name_head = out_dir_COV + "ss_" + number_of_residues;
			path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MAXES_COV.txt";
			Array_IO.write_Double_Array(dist_pca_mode_COV_max, path, 6);
			path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MINS_COV.txt";
			Array_IO.write_Double_Array(dist_pca_mode_COV_min, path, 6);
			path = file_name_head + "_top_" + number_of_modes + "_distance_pca_modes_COV.txt";
			Matrix_IO.write_Matrix(SS_pca_modes, path, 12, 6);
			path = file_name_head + "_top_" + number_of_modes + "_weighted_distance_pca_modes_COV.txt";
			Matrix_IO.write_Matrix(SS_weighted_pca_modes, path, 12, 6);
			path = file_name_head + "_top_" + number_of_modes + "_square_distance_pca_modes_COV.txt";
			Matrix_IO.write_Matrix(SS_pca_square_modes, path, 12, 6);
			path = file_name_head + "_top_" + number_of_modes + "_weighted_square_distance_pca_modes_COV.txt";
			Matrix_IO.write_Matrix(SS_weighted_pca_square_modes, path, 12, 6);
		}

	/* ********************************* CORR METHODS ********************************************************************************** */

	private void get_eigenvalues_CORR()
		{

			double[] ss_evals = evd.getRealEigenvalues();
			trace_CORR = 0.000;
			eigenvalues_CORR = new ArrayList<Double>();
			for (double k : ss_evals)
				{
					eigenvalues_CORR.add(k);
					trace_CORR += k;
				}
			Collections.sort(eigenvalues_CORR, Collections.reverseOrder());
			double max = eigenvalues_CORR.get(0);
			double min = eigenvalues_CORR.get(ROWS_D - 1);
			cond_CORR = Math.abs(max / min);
			file_name_head = out_dir_CORR + "ss_" + number_of_residues;
			List_IO.write_Double_List(eigenvalues_CORR, file_name_head + "_distance_eigenvalues_CORR.txt", 6);
		}

	private void write_top_evals_CORR()
		{
			try
				{
					file_name_head = out_dir_CORR + "ss_" + number_of_residues;
					File top_ss_evals_cov = new File(file_name_head + "_top_" + number_of_modes + "_distance_eigenvalues_CORR.txt");
					BufferedWriter top_ss_evals_writer = new BufferedWriter(new FileWriter(top_ss_evals_cov));
					top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", "Eigenvalue", "% Variance", "Cumulative Variance"));
					top_eigenvalues_CORR = new ArrayList<Double>();
					double cumulative_variance = 0;
					for (int i = 0; i < number_of_modes; i++)
						{
							double val = eigenvalues_CORR.get(i);
							double normed_val = (val / trace_CORR);
							cumulative_variance += normed_val;
							top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", nf.format(val), nf.format(normed_val), nf.format(cumulative_variance)));
							top_eigenvalues_CORR.add(val);
						}
					top_ss_evals_writer.close();

				} catch (IOException io)
				{
					System.err.println("IOException thrown. Could not write the file: " + file_name_head + "_top_" + number_of_modes + "_distance_eigenvalues_CORR.txt");
					io.getMessage();
					io.getStackTrace();
				}
		}

	private void get_top_evects_and_reverse_CORR()
		{

			Matrix ss_evectors = evd.getV();
			Matrix D = evd.getD();
			Matrix precision = ss_evectors.times(D.inverse()).times(ss_evectors.transpose());
			evd = null;
			System.gc();

			top_evectors_dist_CORR = ss_evectors.getMatrix(0, ROWS_D - 1, ROWS_D - number_of_modes, ROWS_D - 1);
			Matrix modes_reversed = new Matrix(ROWS_D, COLS);
			for (int r = 0; r < COLS; r++)
				{
					Matrix col = top_evectors_dist_CORR.getMatrix(0, ROWS_D - 1, COLS - 1 - r, COLS - 1 - r);
					modes_reversed.setMatrix(0, ROWS_D - 1, r, r, col);
				}
			top_evectors_dist_CORR = modes_reversed;

			file_name_head = out_dir_CORR + "ss_" + number_of_residues;
			path = file_name_head + "_top_" + number_of_modes + "_distance_eigenvectors_CORR.txt";
			Matrix_IO.write_Matrix(top_evectors_dist_CORR, path, 12, 6);
			path = file_name_head + "_PRECISION_matrix.txt";
			Matrix_IO.write_Matrix(precision, path, 12, 3);
			ss_evectors = null;
			System.gc();
		}

	private void construct_PCA_Modes_CORR()
		{

			Matrix SS_pca_modes = new Matrix(ROWS_D, number_of_modes);
			Matrix SS_pca_square_modes = new Matrix(ROWS_D, number_of_modes);
			Matrix SS_weighted_pca_modes = new Matrix(ROWS_D, number_of_modes);
			Matrix SS_weighted_pca_square_modes = new Matrix(ROWS_D, number_of_modes);
			dist_pca_mode_CORR_max = new double[number_of_modes];
			dist_pca_mode_CORR_min = new double[number_of_modes];

			for (int a = 0; a < number_of_modes; a++)
				{
					double max = 0;
					double min = 1;

					for (int b = 0; b < ROWS_D; b++)
						{
							double d = top_evectors_dist_CORR.get(b, a);
							double d_abs = Math.abs(d);
							double sq = (d * d);
							double value = top_eigenvalues_CORR.get(a);
							if (value < 0) value = 0;
							double sqrt_val = Math.sqrt(value);
							double wm = sqrt_val * d_abs;
							double w_sq = value * sq;
							SS_pca_modes.set(b, a, d_abs);
							SS_weighted_pca_modes.set(b, a, wm);
							SS_pca_square_modes.set(b, a, sq);
							SS_weighted_pca_square_modes.set(b, a, w_sq);
							if (d_abs >= max)
								{
									max = d_abs;
								}
							if (d_abs <= min)
								{
									min = d_abs;
								}
						}

					dist_pca_mode_CORR_max[a] = max;
					dist_pca_mode_CORR_min[a] = min;
				}
			file_name_head = out_dir_CORR + "ss_" + number_of_residues;
			path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MAXES_CORR.txt";
			Array_IO.write_Double_Array(dist_pca_mode_CORR_max, path, 6);
			path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MINS_CORR.txt";
			Array_IO.write_Double_Array(dist_pca_mode_CORR_min, path, 6);
			path = file_name_head + "_top_" + number_of_modes + "_distance_pca_modes_CORR.txt";
			Matrix_IO.write_Matrix(SS_pca_modes, path, 12, 6);
			path = file_name_head + "_top_" + number_of_modes + "_weighted_distance_pca_modes_CORR.txt";
			Matrix_IO.write_Matrix(SS_weighted_pca_modes, path, 12, 6);
			path = file_name_head + "_top_" + number_of_modes + "_square_distance_pca_modes_CORR.txt";
			Matrix_IO.write_Matrix(SS_pca_square_modes, path, 12, 6);
			path = file_name_head + "_top_" + number_of_modes + "_weighted_square_distance_pca_modes_CORR.txt";
			Matrix_IO.write_Matrix(SS_weighted_pca_square_modes, path, 12, 6);
		}

	/* ********************************** GETTERS ********************************************************************************* */

	public double get_cond_COV()
		{

			return cond_COV;
		}

	public double get_cond_CORR()
		{

			return cond_CORR;
		}

	public Matrix getCOV_dist()
		{

			return COV_dist;
		}

	public Matrix getCORR_dist()
		{

			return CORR_dist;
		}

	public Matrix getResidue_means()
		{

			return residue_means;
		}

	public Matrix getResidue_std_devs()
		{

			return residue_std_devs;
		}

	public double getTrace_COV()
		{

			return trace_COV;
		}

	public double getTrace_CORR()
		{

			return trace_CORR;
		}

	public List<Double> getTop_eigenvalues_COV()
		{

			return top_eigenvalues_COV;
		}

	public List<Double> getTop_eigenvalues_CORR()
		{

			return top_eigenvalues_CORR;
		}

	public Matrix getTop_evectors_dist_COV()
		{

			return top_evectors_dist_COV;
		}

	public Matrix getTop_evectors_dist_CORR()
		{

			return top_evectors_dist_CORR;
		}
}
