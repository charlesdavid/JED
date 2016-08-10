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
 * JED class JED_Get_Distance_Pair_PCA: Gets the COV and CORR PCA for the Residue Distance Pairs subset.
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

public class JED_Get_Distance_Pair_PCA
	{

		String directory, out_dir_dpPCA, out_dir_COV, out_dir_CORR, out_dir_P_CORR, description, file_name_head, path;
		int ROWS_DP, COLS, number_of_modes, number_of_pairs;
		double trace_COV, trace_CORR, trace_P_CORR, cond_COV, cond_CORR, cond_P_CORR, det_COV, det_CORR, det_P_CORR, rank_COV, rank_CORR, rank_P_CORR;
		List<Double> eigenvalues_COV, top_eigenvalues_COV, eigenvalues_CORR, top_eigenvalues_CORR, eigenvalues_P_CORR, top_eigenvalues_P_CORR;
		double[] dist_pca_mode_COV_min, dist_pca_mode_COV_max, dist_pca_mode_CORR_min, dist_pca_mode_CORR_max, dist_pca_mode_P_CORR_min,
				dist_pca_mode_P_CORR_max;
		Matrix distances, centered_data, COV_dist, CORR_dist, pcorr_dist, precision_cov, precision_corr, top_evectors_dist_COV, top_evectors_dist_CORR,
				top_evectors_P_CORR, residue_means, residue_std_devs;
		NumberFormat nf;
		RoundingMode rm;
		EigenvalueDecomposition evd;
		PCA pca;
		boolean success, exist;

		/* ********************************** CONSTRUCTOR *************************************************************************************** */

		/**
		 * Constructor to perform the Distance Pair PCA analysis
		 *
		 * @param dist
		 *            The Distance Pair subset
		 * @param dir
		 *            The working directory
		 * @param des
		 *            The job description
		 * @param modes
		 *            The number of PCA modes to process
		 * @param pairs
		 *            The number of distance pairs
		 */
		JED_Get_Distance_Pair_PCA(Matrix dist, String dir, String des, int modes, int pairs)
			{

				directory = dir;
				description = des;

				out_dir_dpPCA = directory + "JED_RESULTS_" + description + "/dpPCA/";
				exist = new File(out_dir_dpPCA).exists();
				if (!exist) success = (new File(out_dir_dpPCA)).mkdirs();

				out_dir_COV = out_dir_dpPCA + "COV/";
				exist = new File(out_dir_COV).exists();
				if (!exist) success = (new File(out_dir_COV)).mkdirs();

				out_dir_CORR = out_dir_dpPCA + "CORR/";
				exist = new File(out_dir_CORR).exists();
				if (!exist) success = (new File(out_dir_CORR)).mkdirs();

				out_dir_P_CORR = out_dir_dpPCA + "PCORR/";
				exist = new File(out_dir_P_CORR).exists();
				if (!exist) success = (new File(out_dir_P_CORR)).mkdirs();

				distances = dist;
				number_of_pairs = pairs;
				number_of_modes = modes;
				ROWS_DP = number_of_pairs;
				COLS = number_of_modes;

				nf = NumberFormat.getInstance();
				rm = RoundingMode.HALF_UP;
				nf.setRoundingMode(rm);
				nf.setMaximumFractionDigits(3);
				nf.setMinimumFractionDigits(3);
			}

		/* *********************************** DRIVER METHODS *********************************************************************************** */

		/**
		 * Performs the Distance Pair PCA using both COV and CORR (Q and R) models.
		 */
		public void get_Distance_Pair_PCA()
			{

				pca = new PCA(distances);

				distances = null;
				System.gc();

				Do_Cov();
				Do_Corr();
				Do_P_Corr();
			}

		private void Do_Cov()
			{

				COV_dist = pca.get_covariance_matrix_elegant();

				trace_COV = COV_dist.trace();
				cond_COV = COV_dist.cond();
				det_COV = COV_dist.det();
				rank_COV = COV_dist.rank();

				residue_means = pca.getData_means();
				file_name_head = out_dir_dpPCA + number_of_pairs + "_Residue_Pairs";
				path = file_name_head + "_Centroids_of_Variables.txt";
				Matrix_IO.write_Matrix(residue_means, path, 12, 6);

				residue_std_devs = pca.getData_sigmas();
				path = file_name_head + "_Std_Devs_of_Variables.txt";
				Matrix_IO.write_Matrix(residue_std_devs, path, 12, 6);

				file_name_head = out_dir_COV + number_of_pairs + "_Residue_Pairs";
				path = file_name_head + "_COV_matrix.txt";
				Matrix_IO.write_Matrix(COV_dist, path, 12, 6);

				evd = PCA.get_eigenvalue_decomposition(COV_dist);

				get_eigenvalues_COV();
				write_top_evals_COV();
				get_top_evects_and_reverse_COV();
				construct_PCA_Modes_COV();

				System.gc();
			}

		private void Do_Corr()
			{

				CORR_dist = pca.get_R_from_Q(COV_dist);
				COV_dist = null;
				System.gc();

				trace_CORR = CORR_dist.trace();
				cond_CORR = CORR_dist.cond();
				det_CORR = CORR_dist.det();
				rank_CORR = CORR_dist.rank();

				file_name_head = out_dir_CORR + number_of_pairs + "_Residue_Pairs";
				path = file_name_head + "_CORR_matrix.txt";
				Matrix_IO.write_Matrix(CORR_dist, path, 12, 6);

				evd = PCA.get_eigenvalue_decomposition(CORR_dist);

				get_eigenvalues_CORR();
				write_top_evals_CORR();
				get_top_evects_and_reverse_CORR();
				construct_PCA_Modes_CORR();

				System.gc();
			}

		private void Do_P_Corr()
			{

				pcorr_dist = PCA.get_partial_correlation_matrix(precision_cov);

				trace_P_CORR = pcorr_dist.trace();
				cond_P_CORR = pcorr_dist.cond();
				det_P_CORR = pcorr_dist.det();
				rank_P_CORR = pcorr_dist.rank();

				file_name_head = out_dir_P_CORR + "ss_" + number_of_pairs;
				Matrix_IO.write_Matrix(pcorr_dist, file_name_head + "_partial_correlation_matrix.txt", 12, 6);

				evd = PCA.get_eigenvalue_decomposition(pcorr_dist);

				get_eigenvalues_P_CORR();
				write_top_evals_P_CORR();
				get_top_evects_and_reverse_P_CORR();
				construct_PCA_Modes_P_CORR();

				precision_cov = null;
				pcorr_dist = null;
				evd = null;
				System.gc();
			}

		/* ************************************* COV METHODS ************************************************************************************* */

		private void get_eigenvalues_COV()
			{
				double[] ss_evals = evd.getRealEigenvalues();
				eigenvalues_COV = new ArrayList<Double>();
				for (double k : ss_evals)
					{
						eigenvalues_COV.add(k);
					}
				Collections.sort(eigenvalues_COV, Collections.reverseOrder());
				file_name_head = out_dir_COV + number_of_pairs + "_Residue_Pairs";
				path = file_name_head + "_eigenvalues_COV.txt";
				List_IO.write_Double_List(eigenvalues_COV, path, 6);
			}

		private void write_top_evals_COV()
			{
				try
					{
						file_name_head = out_dir_COV + number_of_pairs + "_Residue_Pairs";
						path = file_name_head + "_top_" + number_of_modes + "_eigenvalues_COV.txt";
						File top_ss_evals_cov = new File(path);
						BufferedWriter top_ss_evals_writer = new BufferedWriter(new FileWriter(top_ss_evals_cov));
						top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", "Eigenvalue", "% Variance", "Cumulative Variance"));
						top_eigenvalues_COV = new ArrayList<Double>();
						double cumulative_variance = 0;
						for (int i = 0; i < number_of_modes; i++)
							{
								double val = eigenvalues_COV.get(i);
								double normed_val = (val / trace_COV);
								cumulative_variance += normed_val;
								top_ss_evals_writer
										.write(String.format("%-16s%-16s%-16s%n", nf.format(val), nf.format(normed_val), nf.format(cumulative_variance)));
								top_eigenvalues_COV.add(val);
							}
						top_ss_evals_writer.close();

					} catch (IOException io)
					{
						System.err.println("IOException thrown. Could not write the file: " + path);
						io.printStackTrace();
					}

			}

		private void get_top_evects_and_reverse_COV()
			{

				Matrix ss_evectors = evd.getV();
				Matrix D = evd.getD();
				precision_cov = ss_evectors.times(D.inverse()).times(ss_evectors.transpose());
				evd = null;

				top_evectors_dist_COV = ss_evectors.getMatrix(0, ROWS_DP - 1, ROWS_DP - number_of_modes, ROWS_DP - 1);
				Matrix modes_reversed = new Matrix(ROWS_DP, COLS);
				for (int r = 0; r < COLS; r++)
					{
						Matrix col = top_evectors_dist_COV.getMatrix(0, ROWS_DP - 1, COLS - 1 - r, COLS - 1 - r);
						modes_reversed.setMatrix(0, ROWS_DP - 1, r, r, col);
					}
				top_evectors_dist_COV = modes_reversed;

				file_name_head = out_dir_COV + number_of_pairs + "_Residue_Pairs";
				path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_COV.txt";
				Matrix_IO.write_Matrix(top_evectors_dist_COV, path, 12, 6);
				path = file_name_head + "_inverse_covariance_matrix.txt";
				Matrix_IO.write_Matrix(precision_cov, path, 12, 3);
				ss_evectors = null;
			}

		private void construct_PCA_Modes_COV()
			{

				Matrix SS_pca_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_pca_square_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_weighted_pca_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_weighted_pca_square_modes = new Matrix(ROWS_DP, number_of_modes);
				dist_pca_mode_COV_max = new double[number_of_modes];
				dist_pca_mode_COV_min = new double[number_of_modes];

				for (int a = 0; a < number_of_modes; a++)
					{
						double max = 0;
						double min = 1;

						for (int b = 0; b < ROWS_DP; b++)
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
								if (d_abs >= max) max = d_abs;
								if (d_abs <= min) min = d_abs;
							}

						dist_pca_mode_COV_max[a] = max;
						dist_pca_mode_COV_min[a] = min;
					}
				file_name_head = out_dir_COV + number_of_pairs + "_Residue_Pairs";
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MAXES_COV.txt";
				Array_IO.write_Double_Array(dist_pca_mode_COV_max, path, 6);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MINS_COV.txt";
				Array_IO.write_Double_Array(dist_pca_mode_COV_min, path, 6);
				path = file_name_head + "_top_" + number_of_modes + "_pca_modes_COV.txt";
				Matrix_IO.write_Matrix(SS_pca_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_COV.txt";
				Matrix_IO.write_Matrix(SS_weighted_pca_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_COV.txt";
				Matrix_IO.write_Matrix(SS_pca_square_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_COV.txt";
				Matrix_IO.write_Matrix(SS_weighted_pca_square_modes, path, 12, 6);
			}

		/* ************************************** CORR METHODS *********************************************************************************** */

		private void get_eigenvalues_CORR()
			{
				double[] ss_evals = evd.getRealEigenvalues();
				eigenvalues_CORR = new ArrayList<Double>();
				for (double k : ss_evals)
					{
						eigenvalues_CORR.add(k);
					}
				Collections.sort(eigenvalues_CORR, Collections.reverseOrder());
				file_name_head = out_dir_CORR + number_of_pairs + "_Residue_Pairs";
				path = file_name_head + "_eigenvalues_CORR.txt";
				List_IO.write_Double_List(eigenvalues_CORR, path, 6);
			}

		private void write_top_evals_CORR()
			{
				try
					{
						file_name_head = out_dir_CORR + number_of_pairs + "_Residue_Pairs";
						path = file_name_head + "_top_" + number_of_modes + "_eigenvalues_CORR.txt";
						File top_ss_evals_cov = new File(path);
						BufferedWriter top_ss_evals_writer = new BufferedWriter(new FileWriter(top_ss_evals_cov));
						top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", "Eigenvalue", "% Variance", "Cumulative Variance"));
						top_eigenvalues_CORR = new ArrayList<Double>();
						double cumulative_variance = 0;
						for (int i = 0; i < number_of_modes; i++)
							{
								double val = eigenvalues_CORR.get(i);
								double normed_val = (val / trace_CORR);
								cumulative_variance += normed_val;
								top_ss_evals_writer
										.write(String.format("%-16s%-16s%-16s%n", nf.format(val), nf.format(normed_val), nf.format(cumulative_variance)));
								top_eigenvalues_CORR.add(val);
							}
						top_ss_evals_writer.close();

					} catch (IOException io)
					{
						System.err.println("IOException thrown. Could not read the file: " + path);
						io.printStackTrace();
					}
			}

		private void get_top_evects_and_reverse_CORR()
			{

				Matrix ss_evectors = evd.getV();
				Matrix D = evd.getD();
				precision_corr = ss_evectors.times(D.inverse()).times(ss_evectors.transpose());
				evd = null;
				System.gc();

				top_evectors_dist_CORR = ss_evectors.getMatrix(0, ROWS_DP - 1, ROWS_DP - number_of_modes, ROWS_DP - 1);
				Matrix modes_reversed = new Matrix(ROWS_DP, COLS);
				for (int r = 0; r < COLS; r++)
					{
						Matrix col = top_evectors_dist_CORR.getMatrix(0, ROWS_DP - 1, COLS - 1 - r, COLS - 1 - r);
						modes_reversed.setMatrix(0, ROWS_DP - 1, r, r, col);
					}
				top_evectors_dist_CORR = modes_reversed;

				file_name_head = out_dir_CORR + number_of_pairs + "_Residue_Pairs";
				path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_CORR.txt";
				Matrix_IO.write_Matrix(top_evectors_dist_CORR, path, 12, 6);
				path = file_name_head + "_inverse_correlation_matrix.txt";
				Matrix_IO.write_Matrix(precision_corr, path, 12, 3);
				ss_evectors = null;
				System.gc();
			}

		private void construct_PCA_Modes_CORR()
			{

				Matrix SS_pca_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_pca_square_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_weighted_pca_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_weighted_pca_square_modes = new Matrix(ROWS_DP, number_of_modes);
				dist_pca_mode_CORR_max = new double[number_of_modes];
				dist_pca_mode_CORR_min = new double[number_of_modes];

				for (int a = 0; a < number_of_modes; a++)
					{
						double max = 0;
						double min = 1;

						for (int b = 0; b < ROWS_DP; b++)
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
								if (d_abs >= max) max = d_abs;
								if (d_abs <= min) min = d_abs;
							}

						dist_pca_mode_CORR_max[a] = max;
						dist_pca_mode_CORR_min[a] = min;
					}
				file_name_head = out_dir_CORR + number_of_pairs + "_Residue_Pairs";
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MAXES_CORR.txt";
				Array_IO.write_Double_Array(dist_pca_mode_CORR_max, path, 6);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MINS_CORR.txt";
				Array_IO.write_Double_Array(dist_pca_mode_CORR_min, path, 6);
				path = file_name_head + "_top_" + number_of_modes + "_pca_modes_CORR.txt";
				Matrix_IO.write_Matrix(SS_pca_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_CORR.txt";
				Matrix_IO.write_Matrix(SS_weighted_pca_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_CORR.txt";
				Matrix_IO.write_Matrix(SS_pca_square_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_CORR.txt";
				Matrix_IO.write_Matrix(SS_weighted_pca_square_modes, path, 12, 6);
			}

		/* *********************************** P_CORR METHODS ******************************************************************** */

		private void get_eigenvalues_P_CORR()
			{

				double[] ss_evals = evd.getRealEigenvalues();
				eigenvalues_P_CORR = new ArrayList<Double>();
				for (double k : ss_evals)
					{
						eigenvalues_P_CORR.add(k);
					}
				file_name_head = out_dir_P_CORR + "ss_" + number_of_pairs;
				List_IO.write_Double_List(eigenvalues_P_CORR, file_name_head + "_eigenvalues_P_CORR.txt", 12);
			}

		private void write_top_evals_P_CORR()
			{
				try
					{
						file_name_head = out_dir_P_CORR + "ss_" + number_of_pairs;
						File top_ss_evals_cov = new File(file_name_head + "_top_" + number_of_modes + "_eigenvalues_P_CORR.txt");
						BufferedWriter top_ss_evals_writer = new BufferedWriter(new FileWriter(top_ss_evals_cov));
						top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", "Eigenvalue", "% Variance", "Cumulative Variance"));
						top_eigenvalues_P_CORR = new ArrayList<Double>();
						double cumulative_variance = 0;
						for (int i = 0; i < number_of_modes; i++)
							{
								double val = eigenvalues_P_CORR.get(i);
								double normed_val = (val / trace_P_CORR);
								cumulative_variance += normed_val;
								top_ss_evals_writer
										.write(String.format("%-16s%-16s%-16s%n", nf.format(val), nf.format(normed_val), nf.format(cumulative_variance)));
								top_eigenvalues_P_CORR.add(val);
							}
						top_ss_evals_writer.close();
					} catch (IOException io)
					{
						System.err.println("Could not write to the file: " + file_name_head + "_top_" + number_of_modes + "_eigenvalues_P_CORR.txt");
						io.getMessage();
						io.getStackTrace();
					}
			}

		private void get_top_evects_and_reverse_P_CORR()
			{

				Matrix ss_evectors = evd.getV();

				evd = null;
				System.gc();

				top_evectors_P_CORR = ss_evectors.getMatrix(0, ROWS_DP - 1, 0, number_of_modes - 1);
				file_name_head = out_dir_P_CORR + "ss_" + number_of_pairs;
				path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_P_CORR.txt";
				Matrix_IO.write_Matrix(top_evectors_P_CORR, path, 12, 6);
				ss_evectors = null;
				System.gc();
			}

		private void construct_PCA_Modes_P_CORR()
			{

				Matrix SS_pca_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_pca_square_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_weighted_pca_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_weighted_pca_square_modes = new Matrix(ROWS_DP, number_of_modes);
				dist_pca_mode_P_CORR_max = new double[number_of_modes];
				dist_pca_mode_P_CORR_min = new double[number_of_modes];

				for (int a = 0; a < number_of_modes; a++)
					{
						double max = 0;
						double min = 1;

						for (int b = 0; b < ROWS_DP; b++)
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
								if (d_abs >= max) max = d_abs;
								if (d_abs <= min) min = d_abs;
							}
						dist_pca_mode_P_CORR_max[a] = max;
						dist_pca_mode_P_CORR_min[a] = min;
					}
				file_name_head = out_dir_P_CORR + "ss_" + number_of_pairs;
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MAXES_P_CORR.txt";
				Array_IO.write_Double_Array(dist_pca_mode_P_CORR_max, path, 6);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MINS_P_CORR.txt";
				Array_IO.write_Double_Array(dist_pca_mode_P_CORR_min, path, 6);
				path = file_name_head + "_top_" + number_of_modes + "_pca_modes_P_CORR.txt";
				Matrix_IO.write_Matrix(SS_pca_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_P_CORR.txt";
				Matrix_IO.write_Matrix(SS_weighted_pca_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_P_CORR.txt";
				Matrix_IO.write_Matrix(SS_pca_square_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_P_CORR.txt";
				Matrix_IO.write_Matrix(SS_weighted_pca_square_modes, path, 12, 6);
			}

		/* ************************************** GETTERS ***************************************************************************************** */

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

		public double getTrace_P_CORR()
			{
				return trace_P_CORR;
			}

		public double getDet_COV()
			{
				return det_COV;
			}

		public double getDet_CORR()
			{
				return det_CORR;
			}

		public double getDet_P_CORR()
			{
				return det_P_CORR;
			}

		public double getRank_COV()
			{
				return rank_COV;
			}

		public double getRank_CORR()
			{
				return rank_CORR;
			}

		public double getRank_P_CORR()
			{
				return rank_P_CORR;
			}

		public List<Double> getEigenvalues_COV()
			{
				return eigenvalues_COV;
			}

		public List<Double> getEigenvalues_CORR()
			{
				return eigenvalues_CORR;
			}

		public List<Double> getEigenvalues_P_CORR()
			{
				return eigenvalues_P_CORR;
			}

		public List<Double> getTop_eigenvalues_COV()
			{

				return top_eigenvalues_COV;
			}

		public List<Double> getTop_eigenvalues_CORR()
			{

				return top_eigenvalues_CORR;
			}

		public double[] get_dist_pca_mode_COV_min()
			{

				return dist_pca_mode_COV_min;
			}

		public double[] get_dist_pca_mode_COV_max()
			{

				return dist_pca_mode_COV_max;
			}

		public double[] get_dist_pca_mode_CORR_min()
			{

				return dist_pca_mode_CORR_min;
			}

		public double[] get_dist_pca_mode_CORR_max()
			{

				return dist_pca_mode_CORR_max;
			}

		public Matrix getTop_evectors_dist_COV()
			{

				return top_evectors_dist_COV;
			}

		public Matrix getTop_evectors_dist_CORR()
			{

				return top_evectors_dist_CORR;
			}

		/**
		 * @return the top_eigenvalues_P_CORR
		 */
		public List<Double> getTop_eigenvalues_P_CORR()
			{
				return top_eigenvalues_P_CORR;
			}

		/**
		 * @return the dist_pca_mode_P_CORR_min
		 */
		public double[] getDist_pca_mode_P_CORR_min()
			{
				return dist_pca_mode_P_CORR_min;
			}

		/**
		 * @return the dist_pca_mode_P_CORR_max
		 */
		public double[] getDist_pca_mode_P_CORR_max()
			{
				return dist_pca_mode_P_CORR_max;
			}

		/**
		 * @return the top_evectors_P_CORR
		 */
		public Matrix getTop_evectors_P_CORR()
			{
				return top_evectors_P_CORR;
			}

	}
