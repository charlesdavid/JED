package jed;

import java.io.File;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

import Jama.Matrix;

/**
 * JED class JED_Do_Cartesian: Top class for implementing the Cartesian
 * analysis.
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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * @author Dr. Charles David
 */

public class JED_Do_Cartesian
	{

		public String directory, out_dir, description, pdb_ref_file, rl_SS, Q = "COV", R = "CORR", PC = "P_CORR";
		public int number_of_Chains, number_of_residues, reference_column, number_of_modes_SS;
		public double z_cutoff, percent, trace_COV, trace_CORR, trace_P_CORR, cond_COV, cond_CORR, cond_P_CORR, det_COV, det_CORR, det_P_CORR, rank_COV,
				rank_CORR, rank_P_CORR;
		public List<Integer> residue_list, residue_list_original, residues_read;
		public StringTokenizer sToken;
		public List<String> chain_idents, chain_idents_read, residue_ID_pairs_read;
		public List<Double> transformed_conformation_rmsds, transformed_residue_rmsd_list, top_cartesian_eigenvalues_COV, top_cartesian_eigenvalues_CORR,
				top_cartesian_eigenvalues_P_CORR, eigenvalues_COV, eigenvalues_CORR, eigenvalues_P_CORR;
		public int[] res_list;
		public double[] pca_mode_max_COV, pca_mode_min_COV, pca_mode_max_CORR, pca_mode_min_CORR, pca_mode_max_P_CORR, pca_mode_min_P_CORR;
		public Matrix original_PDB_coordinates, subset_PDB_coordinates, transformed_PDB_coordinates, cov, top_cartesian_evectors_COV, square_pca_modes_COV,
				weighted_square_pca_modes_COV, weighted_pca_modes_COV, normed_projections_COV, projections_COV, corr, top_cartesian_evectors_CORR,
				square_pca_modes_CORR, weighted_square_pca_modes_CORR, weighted_pca_modes_CORR, pca_modes_COV, pca_modes_CORR, normed_projections_CORR,
				projections_CORR, top_cartesian_evectors_P_CORR, square_pca_modes_P_CORR, weighted_square_pca_modes_P_CORR, weighted_pca_modes_P_CORR,
				pca_modes_P_CORR, normed_projections_P_CORR, projections_P_CORR, conf_z_scores, var_z_scores;
		public Vector<Atom> atoms;
		public Matrix trimmmed_PDB_coordinates_COLS, adjusted_PDB_coordinates_rows;

		/* ******************************************** CONSTRUCTORS ********************************************************* */

		/**
		 * Initiates the Cartesian Subset Analysis:
		 *
		 * @param directory
		 *            The working directory
		 * @param description
		 *            The job description
		 * @param pdb_ref
		 *            The PDB Reference File
		 * @param rl_SS
		 *            The Cartesian Residue List
		 * @param reference_column
		 *            The Reference column in the matrix of coordinates.
		 * @param number_of_modes_SS
		 *            The number of Cartesian PCA modes to compute.
		 * @param original_PDB_coordinates
		 *            The Coordinates Matrix
		 */
		@SuppressWarnings("unused")
		public JED_Do_Cartesian(String directory, String description, String pdb_ref, String rl_SS, int reference_column, int number_of_modes_SS,
				Matrix original_PDB_coordinates)
			{

				super();
				this.directory = directory;
				this.description = description;
				out_dir = directory + "JED_RESULTS_" + description + "/cPCA/";
				boolean exist = new File(out_dir).exists();
				if (!exist)
					{
						boolean success = (new File(out_dir)).mkdirs();
					}
				this.pdb_ref_file = pdb_ref;
				this.rl_SS = rl_SS;
				this.reference_column = reference_column;
				this.number_of_modes_SS = number_of_modes_SS;
				this.original_PDB_coordinates = original_PDB_coordinates;
			}

		/* ************************************** DRIVER METHODS ******************************************************* */

		/**
		 * Method for Single Chain PDb files with no chain IDs.
		 */
		public void do_Cartesian()
			{
				read_Residue_List_Single();
				get_Subset();
				get_Transformed_Coords();
				get_Edited_Subset_PDB_Single();
				get_Cartesian_PCA();
				get_Cartesian_DVPs();
			}

		/**
		 * Method for Multi Chain PDb files with chain IDs.
		 */
		public void do_Cartesian_Multi()
			{
				read_Residue_List_Multi();
				get_Subset();
				get_Transformed_Coords();
				get_Edited_Subset_PDB_Multi();
				get_Cartesian_PCA();
				get_Cartesian_DVPs();
			}

		/* ********************************************* METHODS ******************************************************* */

		private void read_Residue_List_Single()
			{
				JED_Read_Residue_List rrl = new JED_Read_Residue_List(directory, description, rl_SS, pdb_ref_file, original_PDB_coordinates);
					{
						rrl.Read_Residue_List_Single();
						residues_read = rrl.getResidues_read();
						residue_list_original = rrl.getResidue_list_original();
						residue_list = rrl.getResidue_list();
						number_of_residues = rrl.getNumber_of_residues();
						res_list = rrl.getRes_list();
					}
				if (number_of_modes_SS > 3 * number_of_residues)
					{
						System.err.println("FATAL ERROR!");
						System.err.println("Number of Cartesian Modes REQUESTED: " + number_of_modes_SS);
						System.err.println("Number of Cartesian Modes AVAILABLE: " + 3 * number_of_residues);
						System.err.println("Possible number of Cartesial Modes is always <= 3*Number_of_Residues in the chosen subset.");
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
			}

		private void read_Residue_List_Multi()
			{
				JED_Read_Residue_List rrl = new JED_Read_Residue_List(directory, description, rl_SS, pdb_ref_file, original_PDB_coordinates);
					{
						rrl.Read_Residue_List_Multi();
						residues_read = rrl.getResidues_read();
						chain_idents_read = rrl.getChain_idents_read();
						residue_list_original = rrl.getResidue_list_original();
						residue_list = rrl.getResidue_list();
						chain_idents = rrl.getChain_idents();
						res_list = rrl.getRes_list();
						number_of_residues = rrl.getNumber_of_residues();
					}
				if (number_of_modes_SS > 3 * number_of_residues)
					{
						System.err.println("FATAL ERROR!");
						System.err.println("Number of Cartesian Modes REQUESTED: " + number_of_modes_SS);
						System.err.println("Number of Cartesian Modes AVAILABLE: " + 3 * number_of_residues);
						System.err.println("Possible number of Cartesial Modes is always <= 3*Number_of_Residues in the chosen subset.");
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
			}

		private void get_Subset()
			{
				JED_Get_Subset ss = new JED_Get_Subset(directory, description, original_PDB_coordinates, res_list);
					{
						ss.set_Output_Directory(out_dir);
						subset_PDB_coordinates = ss.get_Subset_Coords();
					}
			}

		private void get_Transformed_Coords()
			{
				JED_Get_Transformed_Coordinates tf_coords = new JED_Get_Transformed_Coordinates(subset_PDB_coordinates, reference_column, directory,
						description);
					{
						tf_coords.set_Output_Directory(out_dir);
						tf_coords.set_percent_cutoff(percent);
						tf_coords.set_z_cutoff(z_cutoff);
						transformed_PDB_coordinates = tf_coords.get_SS_Transformed_coords();
						transformed_conformation_rmsds = tf_coords.get_SS_Conformation_RMSDs();
						transformed_residue_rmsd_list = tf_coords.get_SS_Residue_RMSDs();
						conf_z_scores = tf_coords.get_conf_Z_scores();
						trimmmed_PDB_coordinates_COLS = tf_coords.get_SS_transformed_coordinates_trimmed_COLS();
						adjusted_PDB_coordinates_rows = tf_coords.get_SS_transformed_coordinates_adjusted_ROWS();

						trimmmed_PDB_coordinates_COLS = null;
						System.gc();
					}
			}

		private void get_Edited_Subset_PDB_Single()
			{
				JED_Edited_Subset_PDB ref_pdb = new JED_Edited_Subset_PDB(directory, pdb_ref_file, transformed_residue_rmsd_list, residue_list_original);
					{
						ref_pdb.set_Output_Directory(out_dir);
						ref_pdb.do_edit_pdb();
						atoms = ref_pdb.getAtoms();
					}
			}

		private void get_Edited_Subset_PDB_Multi()
			{
				JED_Edited_Subset_PDB ref_PDB = new JED_Edited_Subset_PDB(directory, pdb_ref_file, transformed_residue_rmsd_list, chain_idents,
						residue_list_original);
					{
						ref_PDB.set_Output_Directory(out_dir);
						ref_PDB.do_edit_pdb_multi();
						atoms = ref_PDB.getAtoms();
					}
			}

		private void get_Cartesian_PCA()
			{
				JED_Get_Cartesian_PCA c_pca = new JED_Get_Cartesian_PCA(adjusted_PDB_coordinates_rows, number_of_modes_SS, directory, description);
					{
						c_pca.get_Cartesian_PCA();

						// COVARIANCE METHOD
						cond_COV = c_pca.get_cond_COV();
						trace_COV = c_pca.get_trace_COV();
						det_COV = c_pca.get_det_COV();
						rank_COV = c_pca.get_rank_COV();
						eigenvalues_COV = c_pca.getEigenvalues_COV();
						top_cartesian_evectors_COV = c_pca.getTop_evectors_COV();
						top_cartesian_eigenvalues_COV = c_pca.getTop_eigenvalues_COV();
						pca_modes_COV = c_pca.getPca_modes_COV();
						square_pca_modes_COV = c_pca.getSquare_pca_modes_COV();
						weighted_square_pca_modes_COV = c_pca.getWeighted_square_pca_modes_COV();
						weighted_pca_modes_COV = c_pca.getWeighted_pca_modes_COV();
						pca_mode_min_COV = c_pca.get_pca_mode_COV_min();
						pca_mode_max_COV = c_pca.get_pca_mode_COV_max();
						// CORRELATION METHOD
						cond_CORR = c_pca.get_cond_CORR();
						trace_CORR = c_pca.get_trace_CORR();
						det_CORR = c_pca.get_det_CORR();
						rank_CORR = c_pca.get_rank_CORR();
						eigenvalues_CORR = c_pca.getEigenvalues_CORR();
						top_cartesian_evectors_CORR = c_pca.getTop_evectors_CORR();
						top_cartesian_eigenvalues_CORR = c_pca.getEigenvalues_CORR();
						pca_modes_CORR = c_pca.getPca_modes_CORR();
						square_pca_modes_CORR = c_pca.getSquare_pca_modes_CORR();
						weighted_square_pca_modes_CORR = c_pca.getWeighted_square_pca_modes_CORR();
						weighted_pca_modes_CORR = c_pca.getWeighted_pca_modes_CORR();
						pca_mode_min_CORR = c_pca.get_pca_mode_CORR_min();
						pca_mode_max_CORR = c_pca.get_pca_mode_CORR_max();
						// PARTIAL CORRELATION METHOD
						cond_P_CORR = c_pca.get_cond_P_CORR();
						trace_P_CORR = c_pca.get_trace_P_CORR();
						det_P_CORR = c_pca.get_det_P_CORR();
						rank_P_CORR = c_pca.get_rank_P_CORR();
						eigenvalues_P_CORR = c_pca.getEigenvalues_P_CORR();
						top_cartesian_evectors_P_CORR = c_pca.getTop_evectors_P_CORR();
						top_cartesian_eigenvalues_P_CORR = c_pca.getEigenvalues_P_CORR();
						pca_modes_P_CORR = c_pca.getPca_modes_P_CORR();
						square_pca_modes_P_CORR = c_pca.getSquare_pca_modes_P_CORR();
						weighted_square_pca_modes_P_CORR = c_pca.getWeighted_square_pca_modes_P_CORR();
						weighted_pca_modes_P_CORR = c_pca.getWeighted_pca_modes_P_CORR();
						pca_mode_min_P_CORR = c_pca.get_pca_mode_P_CORR_min();
						pca_mode_max_P_CORR = c_pca.get_pca_mode_P_CORR_max();
					}
			}

		private void get_Cartesian_DVPs()
			{
				JED_Get_Cartesian_DVPs pcs_cov = new JED_Get_Cartesian_DVPs(transformed_PDB_coordinates, top_cartesian_evectors_COV, eigenvalues_COV,
						reference_column, directory, description, Q);
					{
						pcs_cov.get_Cartesian_DV_Series();
					}
				JED_Get_Cartesian_DVPs pcs_corr = new JED_Get_Cartesian_DVPs(transformed_PDB_coordinates, top_cartesian_evectors_CORR, eigenvalues_CORR,
						reference_column, directory, description, R);
					{
						pcs_corr.get_Cartesian_DV_Series();
					}
				JED_Get_Cartesian_DVPs pcs_pcorr = new JED_Get_Cartesian_DVPs(transformed_PDB_coordinates, top_cartesian_evectors_P_CORR, eigenvalues_P_CORR,
						reference_column, directory, description, PC);
					{
						pcs_pcorr.get_Cartesian_DV_Series();
					}
				transformed_PDB_coordinates = null;
				System.gc();
			}

		/* ************************************** SETTERS ******************************************************* */

		/**
		 * Sets the Z cutoff for outlier processing
		 *
		 * @param z
		 */
		public void set_z_cutoff(double z)
			{
				z_cutoff = z;
			}

		/**
		 * Sets the Percent of frames to remove for outlier processing
		 *
		 * @param p
		 */
		public void set_percent_cutoff(double p)
			{

				percent = p;
			}

		/* ************************************** GETTERS ******************************************************* */

		public double get_z_cut()
			{

				return z_cutoff;
			}

		public double get_percent_cut()
			{

				return percent;
			}

		public double get_Trace_COV()
			{

				return trace_COV;
			}

		public double get_Trace_CORR()
			{

				return trace_CORR;
			}

		public double get_cond_COV()
			{

				return cond_COV;
			}

		public double get_cond_CORR()
			{

				return cond_CORR;
			}

		public double get_Trace_P_CORR()
			{
				return trace_P_CORR;
			}

		public double get_Cond_COV()
			{
				return cond_COV;
			}

		public double get_Cond_CORR()
			{
				return cond_CORR;
			}

		public double get_Cond_P_CORR()
			{
				return cond_P_CORR;
			}

		public double get_Det_COV()
			{
				return det_COV;
			}

		public double get_Det_CORR()
			{
				return det_CORR;
			}

		public double get_Det_P_CORR()
			{
				return det_P_CORR;
			}

		public double get_Rank_COV()
			{
				return rank_COV;
			}

		public double get_Rank_CORR()
			{
				return rank_CORR;
			}

		public double get_Rank_P_CORR()
			{
				return rank_P_CORR;
			}

		public List<Double> getTop_cartesian_eigenvalues_COV()
			{

				return top_cartesian_eigenvalues_COV;
			}

		public List<Double> getTop_cartesian_eigenvalues_CORR()
			{

				return top_cartesian_eigenvalues_CORR;
			}

		public double[] getPca_mode_max_COV()
			{

				return pca_mode_max_COV;
			}

		public double[] getPca_mode_min_COV()
			{

				return pca_mode_min_COV;
			}

		public double[] getPca_mode_max_CORR()
			{

				return pca_mode_max_CORR;
			}

		public double[] getPca_mode_min_CORR()
			{

				return pca_mode_min_CORR;
			}

		public Matrix getCov()
			{

				return cov;
			}

		public Matrix getTop_cartesian_evectors_COV()
			{

				return top_cartesian_evectors_COV;
			}

		public Matrix getWeighted_pca_modes_COV()
			{

				return weighted_pca_modes_COV;
			}

		public Matrix getSquare_pca_modes_COV()
			{

				return square_pca_modes_COV;
			}

		public Matrix getWeighted_square_pca_modes_COV()
			{

				return weighted_square_pca_modes_COV;
			}

		public Matrix getNormed_projections_COV()
			{

				return normed_projections_COV;
			}

		public Matrix getProjections_COV()
			{

				return projections_COV;
			}

		public Matrix getCorr()
			{

				return corr;
			}

		public Matrix getTop_cartesian_evectors_CORR()
			{

				return top_cartesian_evectors_CORR;
			}

		public Matrix getWeighted_pca_modes_CORR()
			{

				return weighted_pca_modes_CORR;
			}

		public Matrix getSquare_pca_modes_CORR()
			{

				return square_pca_modes_CORR;
			}

		public Matrix getWeighted_square_pca_modes_CORR()
			{

				return weighted_square_pca_modes_CORR;
			}

		public Matrix getNormed_projections_CORR()
			{

				return normed_projections_CORR;
			}

		public Matrix getProjections_CORR()
			{

				return projections_CORR;
			}

		public Vector<Atom> get_atoms()
			{

				return atoms;
			}

		public Matrix getSubset_PDB_coordinates()
			{

				return subset_PDB_coordinates;
			}

		public Matrix getTransformed_PDB_coordinates()
			{

				return transformed_PDB_coordinates;
			}

		public int getNumber_of_residues()
			{

				return number_of_residues;
			}

		public List<Double> getTransformed_conformation_rmsds()
			{

				return transformed_conformation_rmsds;
			}

		public List<Double> getTransformed_residue_rmsd_list()
			{

				return transformed_residue_rmsd_list;
			}

		public Matrix getPca_modes_COV()
			{

				return pca_modes_COV;
			}

		public Matrix getPca_modes_CORR()
			{

				return pca_modes_CORR;
			}

		/**
		 * @return the trace_COV
		 */
		public double getTrace_COV()
			{
				return trace_COV;
			}

		/**
		 * @return the trace_CORR
		 */
		public double getTrace_CORR()
			{
				return trace_CORR;
			}

		/**
		 * @return the trace_P_CORR
		 */
		public double getTrace_P_CORR()
			{
				return trace_P_CORR;
			}

		/**
		 * @return the det_COV
		 */
		public double getDet_COV()
			{
				return det_COV;
			}

		/**
		 * @return the det_CORR
		 */
		public double getDet_CORR()
			{
				return det_CORR;
			}

		/**
		 * @return the det_P_CORR
		 */
		public double getDet_P_CORR()
			{
				return det_P_CORR;
			}

		/**
		 * @return the rank_COV
		 */
		public double getRank_COV()
			{
				return rank_COV;
			}

		/**
		 * @return the rank_CORR
		 */
		public double getRank_CORR()
			{
				return rank_CORR;
			}

		/**
		 * @return the rank_P_CORR
		 */
		public double getRank_P_CORR()
			{
				return rank_P_CORR;
			}

		/**
		 * @return the top_cartesian_eigenvalues_P_CORR
		 */
		public List<Double> getTop_cartesian_eigenvalues_P_CORR()
			{
				return top_cartesian_eigenvalues_P_CORR;
			}

		/**
		 * @return the eigenvalues_COV
		 */
		public List<Double> getEigenvalues_COV()
			{
				return eigenvalues_COV;
			}

		/**
		 * @return the eigenvalues_CORR
		 */
		public List<Double> getEigenvalues_CORR()
			{
				return eigenvalues_CORR;
			}

		/**
		 * @return the eigenvalues_P_CORR
		 */
		public List<Double> getEigenvalues_P_CORR()
			{
				return eigenvalues_P_CORR;
			}

		/**
		 * @return the top_cartesian_evectors_P_CORR
		 */
		public Matrix getTop_cartesian_evectors_P_CORR()
			{
				return top_cartesian_evectors_P_CORR;
			}

		/**
		 * @return the pca_mode_max_P_CORR
		 */
		public double[] getPca_mode_max_P_CORR()
			{
				return pca_mode_max_P_CORR;
			}

		/**
		 * @return the pca_mode_min_P_CORR
		 */
		public double[] getPca_mode_min_P_CORR()
			{
				return pca_mode_min_P_CORR;
			}

		/**
		 * @return the square_pca_modes_P_CORR
		 */
		public Matrix getSquare_pca_modes_P_CORR()
			{
				return square_pca_modes_P_CORR;
			}

		/**
		 * @return the residue_list_original
		 */
		public List<Integer> get_residue_list_original()
			{
				return residue_list_original;
			}
	}
