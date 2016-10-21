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
 * JED class JED_Get_Cartesian_PCA: Gets the PCA using COV and CORR for the Cartesian subset. Copyright (C) 2012 Dr. Charles David
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/license>
 *
 * @author Dr. Charles David
 */
public class JED_Get_Cartesian_PCA
{

	public String directory, out_dir_cPCA, out_dir_COV, out_dir_CORR, out_dir_PCORR, description, file_name_head, path;
	public int number_of_modes, number_of_residues, ROWS, COLS;
	public double trace_COV, trace_CORR, trace_PCORR, cond_COV, cond_CORR, cond_PCORR, rank_COV, rank_CORR, rank_PCORR, det_COV, det_CORR, det_PCORR;
	public List<Double> eigenvalues_COV, top_eigenvalues_COV, eigenvalues_CORR, top_eigenvalues_CORR, eigenvalues_PCORR, top_eigenvalues_PCORR;
	public double[] pca_mode_COV_min, pca_mode_COV_max, pca_mode_CORR_min, pca_mode_CORR_max, pca_mode_PCORR_min, pca_mode_PCORR_max;
	static Matrix input_data, centered_input_data, cov, corr, pcorr, inv_cov, inv_corr, rcov, rcorr, rpcorr, r_inv_cov, r_inv_corr, top_evectors_COV, square_pca_modes_COV,
			weighted_square_pca_modes_COV, weighted_pca_modes_COV, top_evectors_CORR, square_pca_modes_CORR, weighted_square_pca_modes_CORR, weighted_pca_modes_CORR, pca_modes_COV,
			pca_modes_CORR, top_evectors_PCORR, square_pca_modes_PCORR, weighted_square_pca_modes_PCORR, weighted_pca_modes_PCORR, pca_modes_PCORR;
	public EigenvalueDecomposition evd;
	public NumberFormat nf;
	public RoundingMode rm;
	public PCA pca;
	public boolean success, exist;

	/**
	 * Constructor to perform the Cartesian PCA
	 *
	 * @param data
	 *            The Cartesian subset coordinates
	 * @param num_modes
	 *            The number of PCA modes to process
	 * @param dir
	 *            The working directory
	 * @param des
	 *            The job description
	 */
	JED_Get_Cartesian_PCA(Matrix data, int num_modes, String dir, String des)
		{

			nf = NumberFormat.getInstance();
			rm = RoundingMode.HALF_UP;
			nf.setRoundingMode(rm);
			nf.setMaximumFractionDigits(3);
			nf.setMinimumFractionDigits(3);

			directory = dir;
			description = des;

			out_dir_cPCA = directory + "JED_RESULTS_" + description + "/cPCA/";
			exist = new File(out_dir_cPCA).exists();
			if (!exist) success = (new File(out_dir_cPCA)).mkdirs();

			out_dir_COV = out_dir_cPCA + "COV/";
			exist = new File(out_dir_COV).exists();
			if (!exist) success = (new File(out_dir_COV)).mkdirs();

			out_dir_CORR = out_dir_cPCA + "CORR/";
			exist = new File(out_dir_CORR).exists();
			if (!exist) success = (new File(out_dir_CORR)).mkdirs();

			out_dir_PCORR = out_dir_cPCA + "PCORR/";
			exist = new File(out_dir_PCORR).exists();
			if (!exist) success = (new File(out_dir_PCORR)).mkdirs();

			number_of_modes = num_modes;
			input_data = data;
			data = null;
			ROWS = input_data.getRowDimension();
			COLS = number_of_modes;
			number_of_residues = (ROWS / 3);
		}

	/* *************************************** PRIMARY METHODS ****************************************************** */

	/**
	 * Performs the COV, CORR, and PCORR PCA methods
	 */
	public void get_Cartesian_PCA()
		{

			pca = new PCA(input_data);

			input_data = null;
			System.gc();

			Do_Cov_PCA();
			Do_Corr_PCA();
			Do_PCorr_PCA();

		}

	private void Do_Cov_PCA()
		{
			cov = pca.get_covariance_matrix_elegant();
			rcov = PCA.get_reduced_C_matrix(cov);
			trace_COV = cov.trace();
			cond_COV = cov.cond();
			det_COV = cov.det();
			rank_COV = cov.rank();

			Matrix centroids = pca.getData_means();
			Matrix sigmas = pca.getData_sigmas();

			file_name_head = out_dir_cPCA + "ss_" + number_of_residues;

			path = file_name_head + "_means_of_variables.txt";
			Matrix_IO.write_Matrix(centroids, path, 12, 6);
			path = file_name_head + "_std_devs_of_variables.txt";
			Matrix_IO.write_Matrix(sigmas, path, 12, 6);

			file_name_head = out_dir_COV + "ss_" + number_of_residues;

			Matrix_IO.write_Matrix(cov, file_name_head + "_covariance_matrix.txt", 12, 6);
			Matrix_IO.write_Matrix(rcov, file_name_head + "_reduced_covariance_matrix.txt", 12, 6);

			evd = PCA.get_eigenvalue_decomposition(cov);
			get_eigenvalues_COV();
			write_top_evals_COV();
			get_top_evects_and_reverse_COV();
			construct_PCA_Modes_COV();

			evd = null;
			System.gc();
		}

	private void Do_Corr_PCA()
		{

			corr = pca.get_R_from_Q(cov);
			rcorr = PCA.get_reduced_C_matrix(corr);
			trace_CORR = corr.trace();

			cov = null;
			System.gc();

			file_name_head = out_dir_CORR + "ss_" + number_of_residues;

			Matrix_IO.write_Matrix(corr, file_name_head + "_correlation_matrix.txt", 12, 6);
			Matrix_IO.write_Matrix(rcorr, file_name_head + "_reduced_correlation_matrix.txt", 12, 6);

			evd = PCA.get_eigenvalue_decomposition(corr);

			get_eigenvalues_CORR();
			write_top_evals_CORR();
			get_top_evects_and_reverse_CORR();
			construct_PCA_Modes_CORR();

			corr = null;
			evd = null;
			System.gc();
		}

	private void Do_PCorr_PCA()
		{

			pcorr = PCA.get_partial_correlation_matrix(inv_cov);
			rpcorr = PCA.get_reduced_DYN_matrix(pcorr);

			trace_PCORR = pcorr.trace();

			file_name_head = out_dir_PCORR + "ss_" + number_of_residues;

			Matrix_IO.write_Matrix(pcorr, file_name_head + "_partial_correlation_matrix.txt", 12, 6);
			Matrix_IO.write_Matrix(rpcorr, file_name_head + "_reduced_partial_correlation_matrix.txt", 12, 6);

			evd = PCA.get_eigenvalue_decomposition(pcorr);

			get_eigenvalues_PCORR();
			write_top_evals_PCORR();
			get_top_evects_PCORR();
			construct_PCA_Modes_PCORR();

			pcorr = null;
			evd = null;
			System.gc();
		}

	/* **************************************** COV METHODS ************************************************************* */

	private void get_eigenvalues_COV()
		{

			double[] ss_evals = evd.getRealEigenvalues();
			eigenvalues_COV = new ArrayList<>();
			for (double k : ss_evals)
			{
				eigenvalues_COV.add(k);
			}
			Collections.sort(eigenvalues_COV, Collections.reverseOrder());
			List_IO.write_Double_List(eigenvalues_COV, file_name_head + "_eigenvalues_COV.txt", 12);
		}

	private void write_top_evals_COV()
		{
			try
			{
				File top_ss_evals_cov = new File(file_name_head + "_top_" + number_of_modes + "_eigenvalues_COV.txt");
				BufferedWriter top_ss_evals_writer = new BufferedWriter(new FileWriter(top_ss_evals_cov));
				top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", "Eigenvalue", "% Variance", "Cumulative Variance"));
				top_eigenvalues_COV = new ArrayList<>();
				double cumulative_variance = 0;
				for (int e = 0; e < number_of_modes; e++)
				{
					double val = eigenvalues_COV.get(e);
					double normed_val = (val / trace_COV);
					cumulative_variance += normed_val;
					top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", nf.format(val), nf.format(normed_val), nf.format(cumulative_variance)));
					top_eigenvalues_COV.add(val);
				}
				top_ss_evals_writer.close();
			} catch (IOException io)
			{
				System.err.println("Could not write to the file: " + file_name_head + "_top_" + number_of_modes + "_eigenvalues_COV.txt");
				io.getMessage();
				io.getStackTrace();
			}
		}

	private void get_top_evects_and_reverse_COV()
		{

			Matrix ss_evectors = evd.getV();
			Matrix D = evd.getD();

			inv_cov = ss_evectors.times(D.inverse()).times(ss_evectors.inverse());
			r_inv_cov = PCA.get_reduced_DYN_matrix(inv_cov);

			top_evectors_COV = ss_evectors.getMatrix(0, ROWS - 1, ROWS - number_of_modes, ROWS - 1);
			Matrix modes_reversed = new Matrix(ROWS, COLS);
			for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_COV.getMatrix(0, ROWS - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS - 1, r, r, col);
			}
			top_evectors_COV = modes_reversed;

			path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_COV.txt";
			Matrix_IO.write_Matrix(top_evectors_COV, path, 12, 6);

			path = file_name_head + "_inverse_covariance_matrix.txt";
			Matrix_IO.write_Matrix(inv_cov, path, 15, 3);

			path = file_name_head + "_reduced_inverse_covariance_matrix.txt";
			Matrix_IO.write_Matrix(inv_cov, path, 15, 3);

			ss_evectors = null;
			D = null;
			System.gc();
		}

	private void construct_PCA_Modes_COV()
		{

			pca_modes_COV = new Matrix(number_of_residues, number_of_modes);
			weighted_pca_modes_COV = new Matrix(number_of_residues, number_of_modes);
			square_pca_modes_COV = new Matrix(number_of_residues, number_of_modes);
			weighted_square_pca_modes_COV = new Matrix(number_of_residues, number_of_modes);
			pca_mode_COV_max = new double[number_of_modes];
			pca_mode_COV_min = new double[number_of_modes];
			for (int a = 0; a < number_of_modes; a++)
			{
				double max = 0;
				double min = 1;
				for (int b = 0; b < number_of_residues; b++)
				{
					double x = top_evectors_COV.get(b, a);
					double y = top_evectors_COV.get(b + number_of_residues, a);
					double z = top_evectors_COV.get((b + 2 * number_of_residues), a);
					double sq = (x * x + y * y + z * z);
					double sqrt_sq = Math.sqrt(sq);
					double value = eigenvalues_COV.get(a);
					double sqrt_val = Math.sqrt(value);
					double w_sq = sq * value;
					double w_m = sqrt_sq * sqrt_val;
					pca_modes_COV.set(b, a, sqrt_sq);
					weighted_pca_modes_COV.set(b, a, w_m);
					square_pca_modes_COV.set(b, a, sq);
					weighted_square_pca_modes_COV.set(b, a, w_sq);
					if (sq > max) max = sq;
					if (sq < min) min = sq;
				}
				pca_mode_COV_max[a] = max;
				pca_mode_COV_min[a] = min;
			}
			path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MAXES_COV.txt";
			Array_IO.write_Double_Array(pca_mode_COV_max, path, 6);
			path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MINS_COV.txt";
			Array_IO.write_Double_Array(pca_mode_COV_min, path, 6);
			path = file_name_head + "_top_" + number_of_modes + "_pca_modes_COV.txt";
			Matrix_IO.write_Matrix(pca_modes_COV, path, 12, 6);
			path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_COV.txt";
			Matrix_IO.write_Matrix(weighted_pca_modes_COV, path, 12, 6);
			path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_COV.txt";
			Matrix_IO.write_Matrix(square_pca_modes_COV, path, 12, 6);
			path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_COV.txt";
			Matrix_IO.write_Matrix(weighted_square_pca_modes_COV, path, 12, 6);

		}

	/* ***************************************** CORR METHODS ************************************************************** */

	private void get_eigenvalues_CORR()
		{

			double[] ss_evals = evd.getRealEigenvalues();
			eigenvalues_CORR = new ArrayList<>();
			for (double k : ss_evals)
			{
				eigenvalues_CORR.add(k);
			}
			Collections.sort(eigenvalues_CORR, Collections.reverseOrder());
			List_IO.write_Double_List(eigenvalues_CORR, file_name_head + "_eigenvalues_CORR.txt", 12);
		}

	private void write_top_evals_CORR()
		{
			try
			{
				File top_ss_evals_cov = new File(file_name_head + "_top_" + number_of_modes + "_eigenvalues_CORR.txt");
				BufferedWriter top_ss_evals_writer = new BufferedWriter(new FileWriter(top_ss_evals_cov));
				top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", "Eigenvalue", "% Variance", "Cumulative Variance"));
				top_eigenvalues_CORR = new ArrayList<>();
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
				System.err.println("Could not write to the file: " + file_name_head + "_top_" + number_of_modes + "_eigenvalues_CORR.txt");
				io.getMessage();
				io.getStackTrace();
			}
		}

	private void get_top_evects_and_reverse_CORR()
		{

			Matrix ss_evectors = evd.getV();
			Matrix D = evd.getD();

			inv_corr = ss_evectors.times(D.inverse()).times(ss_evectors.inverse());
			r_inv_corr = PCA.get_reduced_DYN_matrix(inv_corr);

			top_evectors_CORR = ss_evectors.getMatrix(0, ROWS - 1, ROWS - number_of_modes, ROWS - 1);
			Matrix modes_reversed = new Matrix(ROWS, COLS);
			for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_CORR.getMatrix(0, ROWS - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS - 1, r, r, col);
			}
			top_evectors_CORR = modes_reversed;

			path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_CORR.txt";
			Matrix_IO.write_Matrix(top_evectors_CORR, path, 12, 6);

			path = file_name_head + "_inverse_correlation_matrix.txt";
			Matrix_IO.write_Matrix(inv_corr, path, 15, 3);

			path = file_name_head + "_reduced_inverse_correlation_matrix.txt";
			Matrix_IO.write_Matrix(r_inv_corr, path, 15, 3);

			ss_evectors = null;
			D = null;
			System.gc();
		}

	private void construct_PCA_Modes_CORR()
		{

			pca_modes_CORR = new Matrix(number_of_residues, number_of_modes);
			weighted_pca_modes_CORR = new Matrix(number_of_residues, number_of_modes);
			square_pca_modes_CORR = new Matrix(number_of_residues, number_of_modes);
			weighted_square_pca_modes_CORR = new Matrix(number_of_residues, number_of_modes);
			pca_mode_CORR_max = new double[number_of_modes];
			pca_mode_CORR_min = new double[number_of_modes];
			for (int a = 0; a < number_of_modes; a++)
			{
				double max = 0;
				double min = 1;
				for (int b = 0; b < number_of_residues; b++)
				{
					double x = top_evectors_CORR.get(b, a);
					double y = top_evectors_CORR.get(b + number_of_residues, a);
					double z = top_evectors_CORR.get((b + 2 * number_of_residues), a);
					double sq = (x * x + y * y + z * z);
					double sqrt_sq = Math.sqrt(sq);
					double value = eigenvalues_CORR.get(a);
					double sqrt_val = Math.sqrt(value);
					double w_sq = sq * value;
					double w_m = sqrt_sq * sqrt_val;
					pca_modes_CORR.set(b, a, sqrt_sq);
					weighted_pca_modes_CORR.set(b, a, w_m);
					square_pca_modes_CORR.set(b, a, sq);
					weighted_square_pca_modes_CORR.set(b, a, w_sq);
					if (sq > max) max = sq;
					if (sq < min) min = sq;
				}
				pca_mode_CORR_max[a] = max;
				pca_mode_CORR_min[a] = min;
			}
			path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MAXES_CORR.txt";
			Array_IO.write_Double_Array(pca_mode_CORR_max, path, 6);
			path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MINS_CORR.txt";
			Array_IO.write_Double_Array(pca_mode_CORR_min, path, 6);
			path = file_name_head + "_top_" + number_of_modes + "_pca_modes_CORR.txt";
			Matrix_IO.write_Matrix(pca_modes_CORR, path, 12, 6);
			path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_CORR.txt";
			Matrix_IO.write_Matrix(weighted_pca_modes_CORR, path, 12, 6);
			path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_CORR.txt";
			Matrix_IO.write_Matrix(square_pca_modes_CORR, path, 12, 6);
			path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_CORR.txt";
			Matrix_IO.write_Matrix(weighted_square_pca_modes_CORR, path, 12, 6);
		}

	/* ********************************************* PCORR METHODS ************************************************************** */

	private void get_eigenvalues_PCORR()
		{

			/* Note that the pcorr matrix has -1's along its main diagonal, making it a 'negative definite' matrix */
			/* Therefore, the eigenvalues with the largest absolute value are the smallest (most negative). */

			double[] ss_evals = evd.getRealEigenvalues();
			eigenvalues_PCORR = new ArrayList<>();
			for (double k : ss_evals)
			{
				eigenvalues_PCORR.add(k);
			}
			Collections.sort(eigenvalues_PCORR, Collections.reverseOrder());
			List_IO.write_Double_List(eigenvalues_PCORR, file_name_head + "_eigenvalues_PCORR.txt", 12);
		}

	private void write_top_evals_PCORR()
		{
			try
			{
				File top_ss_evals_cov = new File(file_name_head + "_top_" + number_of_modes + "_eigenvalues_PCORR.txt");
				BufferedWriter top_ss_evals_writer = new BufferedWriter(new FileWriter(top_ss_evals_cov));
				top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", "Eigenvalue", "% Variance", "Cumulative Variance"));
				top_eigenvalues_PCORR = new ArrayList<>();
				double cumulative_variance = 0;
				for (int i = 0; i < number_of_modes; i++)
				{
					double val = eigenvalues_PCORR.get(i);
					double normed_val = (val / trace_PCORR);
					cumulative_variance += normed_val;
					top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", nf.format(val), nf.format(normed_val), nf.format(cumulative_variance)));
					top_eigenvalues_PCORR.add(val);
				}
				top_ss_evals_writer.close();
			} catch (IOException io)
			{
				System.err.println("Could not write to the file: " + file_name_head + "_top_" + number_of_modes + "_eigenvalues_PCORR.txt");
				io.getMessage();
				io.getStackTrace();
			}
		}

	private void get_top_evects_PCORR()
		{
			Matrix ss_evectors = evd.getV();
			evd = null;
			System.gc();

			top_evectors_PCORR = ss_evectors.getMatrix(0, ROWS - 1, ROWS - number_of_modes, ROWS - 1);
			Matrix modes_reversed = new Matrix(ROWS, COLS);
			for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_PCORR.getMatrix(0, ROWS - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS - 1, r, r, col);
			}
			top_evectors_PCORR = modes_reversed;

			path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_PCORR.txt";
			Matrix_IO.write_Matrix(top_evectors_PCORR, path, 12, 6);
			ss_evectors = null;
			System.gc();
		}

	private void construct_PCA_Modes_PCORR()
		{

			pca_modes_PCORR = new Matrix(number_of_residues, number_of_modes);
			weighted_pca_modes_PCORR = new Matrix(number_of_residues, number_of_modes);
			square_pca_modes_PCORR = new Matrix(number_of_residues, number_of_modes);
			weighted_square_pca_modes_PCORR = new Matrix(number_of_residues, number_of_modes);
			pca_mode_PCORR_max = new double[number_of_modes];
			pca_mode_PCORR_min = new double[number_of_modes];
			for (int a = 0; a < number_of_modes; a++)
			{
				double max = 0;
				double min = 1;
				for (int b = 0; b < number_of_residues; b++)
				{
					double x = top_evectors_PCORR.get(b, a);
					double y = top_evectors_PCORR.get(b + number_of_residues, a);
					double z = top_evectors_PCORR.get((b + 2 * number_of_residues), a);
					double sq = (x * x + y * y + z * z);
					double sqrt_sq = Math.sqrt(sq);
					/* Note that the 'top' eigenvalues for PCORR will be close to zero */
					double value = eigenvalues_PCORR.get(a);
					double value_abs = Math.abs(value);
					double sqrt_val = Math.sqrt(value_abs);
					double w_sq = sq * value_abs;
					double w_m = sqrt_sq * sqrt_val;
					pca_modes_PCORR.set(b, a, sqrt_sq);
					weighted_pca_modes_PCORR.set(b, a, w_m);
					square_pca_modes_PCORR.set(b, a, sq);
					weighted_square_pca_modes_PCORR.set(b, a, w_sq);
					if (sq > max) max = sq;
					if (sq < min) min = sq;
				}
				pca_mode_PCORR_max[a] = max;
				pca_mode_PCORR_min[a] = min;
			}
			path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MAXES_PCORR.txt";
			Array_IO.write_Double_Array(pca_mode_PCORR_max, path, 6);
			path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MINS_PCORR.txt";
			Array_IO.write_Double_Array(pca_mode_PCORR_min, path, 6);
			path = file_name_head + "_top_" + number_of_modes + "_pca_modes_PCORR.txt";
			Matrix_IO.write_Matrix(pca_modes_PCORR, path, 12, 6);
			path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_PCORR.txt";
			Matrix_IO.write_Matrix(weighted_pca_modes_PCORR, path, 18, 12);
			path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_PCORR.txt";
			Matrix_IO.write_Matrix(square_pca_modes_PCORR, path, 18, 12);
			path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_PCORR.txt";
			Matrix_IO.write_Matrix(weighted_square_pca_modes_PCORR, path, 18, 12);
		}

	/* ******************************************* GETTERS ***************************************************** */

	/**
	 * Returns the condition number of the covariance matrix
	 *
	 * @return
	 */
	public double get_cond_COV()
		{

			return cond_COV;
		}

	/**
	 * Returns the condition number of the correlation matrix
	 *
	 * @return
	 */
	public double get_cond_CORR()
		{

			return cond_CORR;
		}

	/**
	 * Returns the condition number of the partial correlation matrix
	 *
	 * @return
	 */
	public double get_cond_PCORR()
		{

			return cond_PCORR;
		}

	/**
	 * Returns the Trace of the covariance matrix
	 *
	 * @return
	 */
	public double get_trace_COV()
		{

			return trace_COV;
		}

	/**
	 * Returns the Trace of the correlation matrix
	 *
	 * @return
	 */
	public double get_trace_CORR()
		{

			return trace_CORR;
		}

	/**
	 * Returns the Trace of the partial correlation matrix
	 *
	 * @return
	 */
	public double get_trace_PCORR()
		{

			return trace_PCORR;
		}

	/**
	 * Returns the Eigenvalues of the covariance matrix
	 *
	 * @return
	 */
	public List<Double> getEigenvalues_COV()
		{

			return eigenvalues_COV;
		}

	/**
	 * Returns the Eigenvalues of the correlation matrix
	 *
	 * @return
	 */
	public List<Double> getEigenvalues_CORR()
		{

			return eigenvalues_CORR;
		}

	/**
	 * Returns the Eigenvalues of the partial correlation matrix
	 *
	 * @return
	 */
	public List<Double> getEigenvalues_PCORR()
		{

			return eigenvalues_PCORR;
		}

	/**
	 * Returns the determinant of the covariance matrix
	 *
	 * @return
	 */
	public double get_det_COV()
		{

			return det_COV;
		}

	/**
	 * Returns the determinant of the correlation matrix
	 *
	 * @return
	 */
	public double get_det_CORR()
		{

			return det_CORR;
		}

	/**
	 * Returns the determinant of the partial correlation matrix
	 *
	 * @return
	 */
	public double get_det_PCORR()
		{

			return det_PCORR;
		}

	/**
	 * Returns the rank of the covariance matrix
	 *
	 * @return
	 */
	public double get_rank_COV()
		{

			return rank_COV;
		}

	/**
	 * Returns the rank of the correlation matrix
	 *
	 * @return
	 */
	public double get_rank_CORR()
		{

			return rank_CORR;
		}

	/**
	 * Returns the rank of the partial correlation matrix
	 *
	 * @return
	 */
	public double get_rank_PCORR()
		{

			return rank_PCORR;
		}

	/**
	 * Returns the array of the PCA mode minimums from the COV analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_COV_min()
		{

			return pca_mode_COV_min;
		}

	/**
	 * Returns the array of the PCA mode maximums from the COV analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_COV_max()
		{

			return pca_mode_COV_max;
		}

	/**
	 * Returns the array of the PCA mode minimums from the CORR analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_CORR_min()
		{

			return pca_mode_CORR_min;
		}

	/**
	 * Returns the array of the PCA mode maximums from the CORR analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_CORR_max()
		{

			return pca_mode_CORR_max;
		}

	/**
	 *
	 * /** Returns the array of the PCA mode minimums from the PCORR analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_PCORR_min()
		{

			return pca_mode_PCORR_min;
		}

	/**
	 * Returns the array of the PCA mode maximums from the PCORR analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_PCORR_max()
		{

			return pca_mode_PCORR_max;
		}

	/**
	 * Returns the top eigenvectors from the COV analysis
	 *
	 * @return
	 */
	public Matrix getTop_evectors_COV()
		{

			return top_evectors_COV;
		}

	/**
	 * Returns the Square PCA modes from the COV analysis
	 *
	 * @return
	 */
	public Matrix getSquare_pca_modes_COV()
		{

			return square_pca_modes_COV;
		}

	/**
	 * Returns the Weighted Square PCA modes from the COV analysis
	 *
	 * @return
	 */
	public Matrix getWeighted_square_pca_modes_COV()
		{

			return weighted_square_pca_modes_COV;
		}

	/**
	 * Returns the Weighted PCA modes from the COV analysis
	 *
	 * @return
	 */
	public Matrix getWeighted_pca_modes_COV()
		{

			return weighted_pca_modes_COV;
		}

	/**
	 * Returns the Top eigenvectors from the CORR analysis
	 *
	 * @return
	 */
	public Matrix getTop_evectors_CORR()
		{

			return top_evectors_CORR;
		}

	/**
	 * Returns the Square PCA modes from the CORR analysis
	 *
	 * @return
	 */
	public Matrix getSquare_pca_modes_CORR()
		{

			return square_pca_modes_CORR;
		}

	/**
	 * Returns the Weighted Square PCA modes from the CORR analysis
	 *
	 * @return
	 */
	public Matrix getWeighted_square_pca_modes_CORR()
		{

			return weighted_square_pca_modes_CORR;
		}

	/**
	 * Returns the Weighted PCA modes from the CORR analysis
	 *
	 * @return
	 */
	public Matrix getWeighted_pca_modes_CORR()
		{

			return weighted_pca_modes_CORR;
		}

	/**
	 * Returns the PCA modes from the COV analysis
	 *
	 * @return
	 */
	public Matrix getPca_modes_COV()
		{

			return pca_modes_COV;
		}

	/**
	 * Returns the PCA modes from the CORR analysis
	 *
	 * @return
	 */
	public Matrix getPca_modes_CORR()
		{

			return pca_modes_CORR;
		}

	/**
	 * Returns the Top eigenvalues from the COV analysis
	 *
	 * @return
	 */
	public List<Double> getTop_eigenvalues_COV()
		{

			return top_eigenvalues_COV;
		}

	/**
	 * Returns the Top eigenvalues from the CORR analysis
	 *
	 * @return
	 */
	public List<Double> getTop_eigenvalues_CORR()
		{

			return top_eigenvalues_CORR;
		}

	/**
	 * Returns the Top eigenvectors from the PCORR analysis
	 *
	 * @return
	 */
	public Matrix getTop_evectors_PCORR()
		{

			return top_evectors_PCORR;
		}

	/**
	 * Returns the Square PCA modes from the PCORR analysis
	 *
	 * @return
	 */
	public Matrix getSquare_pca_modes_PCORR()
		{

			return square_pca_modes_PCORR;
		}

	/**
	 * Returns the Weighted Square PCA modes from the PCORR analysis
	 *
	 * @return
	 */
	public Matrix getWeighted_square_pca_modes_PCORR()
		{

			return weighted_square_pca_modes_PCORR;
		}

	/**
	 * Returns the Weighted PCA modes from the PCORR analysis
	 *
	 * @return
	 */
	public Matrix getWeighted_pca_modes_PCORR()
		{

			return weighted_pca_modes_PCORR;
		}

	/**
	 * Returns the PCA modes from the PCORR analysis
	 *
	 * @return
	 */
	public Matrix getPca_modes_PCORR()
		{

			return pca_modes_PCORR;
		}
}
