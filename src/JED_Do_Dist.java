package jed;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

import Jama.Matrix;

/**
 * JED class JED _Do_Dist. Top class for the All-to-All Distance analysis.
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
 * 
 */
public class JED_Do_Dist
{

	public int number_of_residues_dist, reference_column, number_of_modes_dist;
	public String directory, out_dir, description, pdb_ref_file, rl_dist;
	public double trace_dist_COV, trace_dist_CORR, cond_cov, cond_corr, z_cutoff;
	public List<Integer> residue_list_dist, residue_list_dist_orig, residues_read;
	public List<String> chain_idents, chain_idents_read;
	public List<Double> transformed_conformation_rmsds_dist, transformed_residue_rmsd_list_dist, top_distance_eigenvalues_COV, top_distance_eigenvalues_CORR;
	public int[] res_list_dist, res_list_dist_orig;
	public double[] residue_distance_means, residue_distance_std_devs, dist_pca_mode_COV_min, dist_pca_mode_COV_max, dist_pca_mode_CORR_min, dist_pca_mode_CORR_max;
	public StringTokenizer sToken;
	public Vector<Atom> atoms;
	public static Matrix distance_matrix, distance_matrix_cond, original_PDB_coordinates, subset_PDB_coordinates_dist, transformed_subset_PDB_coordinates_dist, cov_dist, corr_dist,
			top_distance_evectors_COV, top_distance_evectors_CORR, normed_projections_dist_COV, projections_dist_COV, projections_dist_CORR, normed_projections_dist_CORR,
			weighted_normed_projections_dist_COV, weighted_normed_projections_dist_CORR, Z_scores, counts;
	NumberFormat nf;
	RoundingMode rm;
	boolean success, exist;

	/* ************************************ CONSTRUCTORS ************************************************************************* */

	/**
	 * Initiates the ALL to ALL Distance Analysis:
	 * 
	 * @param dir
	 *            The working directory
	 * @param desc
	 *            The job description
	 * @param ref_file
	 *            The PDB Reference file
	 * @param res_d
	 *            The list of Residue Pairs to analyze
	 * @param ref_col
	 *            The Reference Column in the Coordinates Matrix
	 * @param modes_dist
	 *            The number of Distance Pair PCA modes to process
	 * @param PDB_coordinates
	 *            The Matrix of Coordinates
	 */
	public JED_Do_Dist(String dir, String desc, String ref_file, String res_d, int ref_col, int modes_dist, Matrix PDB_coordinates)
		{

			super();
			directory = dir;
			description = desc;
			out_dir = directory + "JED_RESULTS_" + description + "/dPCA/";
			exist = new File(out_dir).exists();
			if (!exist)
				{
					success = (new File(out_dir)).mkdirs();
				}
			pdb_ref_file = ref_file;
			rl_dist = res_d;
			reference_column = ref_col;
			number_of_modes_dist = modes_dist;
			original_PDB_coordinates = PDB_coordinates;

			nf = NumberFormat.getInstance();
			rm = RoundingMode.HALF_UP;
			nf.setRoundingMode(rm);
			nf.setMaximumFractionDigits(3);
			nf.setMinimumFractionDigits(3);
		}

	/* **************************************** DRIVER METHODS ********************************************************************* */

	/**
	 * Method for processing Single Chain PDBs with no chain IDs.
	 */
	public void do_Dist()
		{

			read_Residue_List_Single();
			get_Subset();
			get_Subset_PDB_Single();
			get_All_to_All_Distances();
			remove_Outliers_Z_Score();
			get_Distance_PCA();
			get_Distance_DVPs();
			write_distance_stats();
		}

	/**
	 * Method for processing Multi chain PDBs with chain IDs.
	 */
	public void do_Dist_Multi()
		{

			read_Residue_List_Multi();
			get_Subset();
			get_Subset_PDB_Multi();
			get_All_to_All_Distances();
			remove_Outliers_Z_Score();
			get_Distance_PCA();
			get_Distance_DVPs();
			write_distance_stats_multi();
		}

	/* ******************************************* METHODS ************************************************************************ */

	private void get_Distance_DVPs()
		{
			JED_Get_Distance_DVPs SS_dvp_COV = new JED_Get_Distance_DVPs(distance_matrix, top_distance_evectors_COV, top_distance_eigenvalues_COV, reference_column, directory, description, "COV",
					number_of_modes_dist, number_of_residues_dist);
				{
					// COVARIANCE METHOD
					SS_dvp_COV.get_Distance_DV_Series();
					SS_dvp_COV.get_DVPs();
				}
			JED_Get_Distance_DVPs SS_dvp_CORR = new JED_Get_Distance_DVPs(distance_matrix, top_distance_evectors_CORR, top_distance_eigenvalues_CORR, reference_column, directory, description, "CORR",
					number_of_modes_dist, number_of_residues_dist);
				{
					// CORRELATION METHOD
					SS_dvp_CORR.get_Distance_DV_Series();
					SS_dvp_CORR.get_DVPs();
				}
		}

	private void get_Distance_PCA()
		{
			JED_Get_Distance_PCA dpca = new JED_Get_Distance_PCA(distance_matrix_cond, directory, description, number_of_modes_dist, number_of_residues_dist);
				{
					dpca.get_Distance_PCA();

					residue_distance_means = dpca.getResidue_means().getColumnPackedCopy();
					residue_distance_std_devs = dpca.getResidue_std_devs().getColumnPackedCopy();

					// COVARIANCE METHOD
					cond_cov = dpca.get_cond_COV();
					top_distance_evectors_COV = dpca.getTop_evectors_dist_COV();
					top_distance_eigenvalues_COV = dpca.getTop_eigenvalues_COV();
					trace_dist_COV = dpca.getTrace_COV();
					dist_pca_mode_COV_max = dpca.dist_pca_mode_COV_max;
					dist_pca_mode_COV_min = dpca.dist_pca_mode_COV_min;

					// CORRELATION METHOD
					cond_corr = dpca.get_cond_CORR();
					top_distance_evectors_CORR = dpca.getTop_evectors_dist_CORR();
					top_distance_eigenvalues_CORR = dpca.getTop_eigenvalues_CORR();
					trace_dist_CORR = dpca.getTrace_CORR();
					dist_pca_mode_CORR_max = dpca.dist_pca_mode_CORR_max;
					dist_pca_mode_CORR_min = dpca.dist_pca_mode_CORR_min;
				}
		}

	private void remove_Outliers_Z_Score()
		{
			if (z_cutoff > 0)
				{
					Adjust_Outliers_by_Z_Score cv = new Adjust_Outliers_by_Z_Score(distance_matrix);
						{
							cv.set_Z_threshold(z_cutoff);
							cv.adjust_row_data();
							distance_matrix_cond = cv.get_coorinates_adjusted();
							Z_scores = cv.get_z_scores();
							Matrix_IO.write_Matrix(Z_scores, out_dir + "ss_" + number_of_residues_dist + "_distance_Z_scores.txt", 6, 1);
							counts = cv.get_counts();
							Matrix_IO.write_Matrix(counts, out_dir + "ss_" + number_of_residues_dist + "_outliers_per_variable.txt", 6, 0);
							// residue_distance_means = cv.get_means();
							// residue_distance_std_devs = cv.get_std_deviations();
						}
				}
			if (z_cutoff == 0) distance_matrix_cond = distance_matrix;
		}

	private void get_All_to_All_Distances()
		{
			JED_Get_All_to_All_Distances gaad = new JED_Get_All_to_All_Distances(subset_PDB_coordinates_dist, directory, description);
				{
					distance_matrix = gaad.get_All_To_All_Distances();
				}
		}

	private void get_Subset_PDB_Single()
		{
			JED_Edited_Subset_PDB ref_pdb = new JED_Edited_Subset_PDB(residue_list_dist_orig, directory, pdb_ref_file);
				{
					ref_pdb.set_Output_Directory(out_dir);
					ref_pdb.get_pdb_subset();
					atoms = ref_pdb.getAtoms();
				}
		}

	private void get_Subset()
		{
			JED_Get_Subset dss = new JED_Get_Subset(directory, description, original_PDB_coordinates, res_list_dist);
				{
					dss.set_Output_Directory(out_dir);
					subset_PDB_coordinates_dist = dss.get_Subset_Coords();
				}
		}

	private void read_Residue_List_Single()
		{
			JED_Read_Residue_List rrl = new JED_Read_Residue_List(directory, description, rl_dist, pdb_ref_file, original_PDB_coordinates);
				{
					rrl.Read_Residue_List_Single();
					residues_read = rrl.getResidues_read();
					residue_list_dist = rrl.getResidue_list();
					residue_list_dist_orig = rrl.getResidue_list_original();
					number_of_residues_dist = residue_list_dist.size();
					res_list_dist = rrl.getRes_list();
					res_list_dist_orig = rrl.getRes_list_orig();
				}
			if (number_of_modes_dist > (number_of_residues_dist * (number_of_residues_dist - 1) / 2))
				{
					System.err.println("FATAL ERROR!");
					System.err.println("Number of Cartesian Modes REQUESTED: " + number_of_modes_dist);
					System.err.println("Number of Cartesian Modes AVAILABLE: " + (number_of_residues_dist * (number_of_residues_dist - 1) / 2));
					System.err.println("Possible number of Cartesial Modes is always <= (Number_of_Residues_Dist)(Number_of_Residues_Dist-1)/2.");
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
		}

	private void write_distance_stats()
		{
			try
				{
					File d_stats = new File(out_dir + "ss_" + res_list_dist.length + "_distance_residues_stats.txt");
					BufferedWriter d_stats_writer = new BufferedWriter(new FileWriter(d_stats));
					d_stats_writer.write("MEANs and STANDARD DEVIATIONs for the Residue Distances: " + "\n");
					d_stats_writer.write(String.format("%-12s%-16s%-16s%-16s%n", "Res1", "Res2", "Mean", "Std_Dev"));
					int j, k = 0;
					for (int i = 0; i < number_of_residues_dist; i++)
						{
							for (j = (i + 1); j < number_of_residues_dist; j++)
								{
									d_stats_writer.write(String.format("%-12s%-16s%-16s%-16s%n", res_list_dist_orig[i], res_list_dist_orig[j], nf.format(residue_distance_means[k]),
											nf.format(residue_distance_std_devs[k])));
									k++;
								}
						}
					d_stats_writer.close();

				} catch (IOException io)
				{
					System.err.println("IOException thrown. Could not write the file: " + out_dir + "ss_" + res_list_dist.length + "_distance_residues_stats.txt");
					System.err.println("Terminating program execution.");
					io.printStackTrace();
					System.exit(0);
				}
		}

	private void get_Subset_PDB_Multi()
		{
			JED_Edited_Subset_PDB ref_pdb = new JED_Edited_Subset_PDB(chain_idents, residue_list_dist_orig, directory, pdb_ref_file);
				{
					ref_pdb.set_Output_Directory(out_dir);
					ref_pdb.get_pdb_subset_multi();
					atoms = ref_pdb.getAtoms();
				}
		}

	private void read_Residue_List_Multi()
		{
			JED_Read_Residue_List rrl = new JED_Read_Residue_List(directory, description, rl_dist, pdb_ref_file, original_PDB_coordinates);
				{
					rrl.Read_Residue_List_Multi();
					residues_read = rrl.getResidues_read();
					chain_idents_read = rrl.getChain_idents_read();
					chain_idents = rrl.getChain_idents();
					residue_list_dist = rrl.getResidue_list();
					residue_list_dist_orig = rrl.getResidue_list_original();
					number_of_residues_dist = residue_list_dist.size();
					res_list_dist = rrl.getRes_list();
					res_list_dist_orig = rrl.getRes_list_orig();
				}
			if (number_of_modes_dist > (number_of_residues_dist * (number_of_residues_dist - 1) / 2))
				{
					System.err.println("FATAL ERROR!");
					System.err.println("Number of Cartesian Modes REQUESTED: " + number_of_modes_dist);
					System.err.println("Number of Cartesian Modes AVAILABLE: " + (number_of_residues_dist * (number_of_residues_dist - 1) / 2));
					System.err.println("Possible number of Cartesial Modes is always <= (Number_of_Residues_Dist)(Number_of_Residues_Dist-1)/2.");
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
		}

	private void write_distance_stats_multi()
		{
			try
				{
					File d_stats = new File(out_dir + "ss_" + res_list_dist.length + "_distance_residues_stats.txt");
					BufferedWriter d_stats_writer = new BufferedWriter(new FileWriter(d_stats));
					d_stats_writer.write("MEANs and STANDARD DEVIATIONs for the Residue Distances: " + "\n");
					d_stats_writer.write(String.format("%-12s%-16s%-16s%-16s%n", "Res1", "Res2", "Mean", "Std_Dev"));
					int j, k = 0;
					for (int i = 0; i < number_of_residues_dist; i++)
						{
							for (j = (i + 1); j < number_of_residues_dist; j++)
								{
									d_stats_writer.write(String.format("%-12s%-16s%-16s%-16s%n", chain_idents.get(i) + res_list_dist_orig[i], chain_idents.get(j) + res_list_dist_orig[j],
											nf.format(residue_distance_means[k]), nf.format(residue_distance_std_devs[k])));
									k++;
								}
						}
					d_stats_writer.close();

				} catch (IOException io)
				{
					System.err.println("IOException thrown. Could not write the file: " + out_dir + "ss_" + res_list_dist.length + "_distance_residues_stats.txt");
					System.err.println("Terminating program execution.");
					io.printStackTrace();
					System.exit(0);
				}
		}

	/* ************************************* SETTERS *************************************************************************** */

	/**
	 * Sets the Z cutoff for outlier processing.
	 * 
	 * @param z
	 *            The Z cutoff
	 */
	public void set_z_cutoff(double z)
		{

			z_cutoff = z;
		}

	/* ************************************* GETTERS *************************************************************************** */

	public double getTrace_dist_COV()
		{

			return trace_dist_COV;
		}

	public double getTrace_dist_CORR()
		{

			return trace_dist_CORR;
		}

	public double get_cond_cov()
		{

			return cond_cov;
		}

	public double get_cond_corr()
		{

			return cond_corr;
		}

	public List<Integer> getResidue_list_dist()
		{

			return residue_list_dist;
		}

	public List<Integer> getResidue_list_dist_orig()
		{

			return residue_list_dist_orig;
		}

	public List<Double> getTransformed_conformation_rmsds_dist()
		{

			return transformed_conformation_rmsds_dist;
		}

	public List<Double> getTransformed_residue_rmsd_list_dist()
		{

			return transformed_residue_rmsd_list_dist;
		}

	public List<Double> getTop_distance_eigenvalues_COV()
		{

			return top_distance_eigenvalues_COV;
		}

	public List<Double> getTop_distance_eigenvalues_CORR()
		{

			return top_distance_eigenvalues_CORR;
		}

	public double[] getResidue_distance_means()
		{

			return residue_distance_means;
		}

	public double[] getResidue_distance_std_devs()
		{

			return residue_distance_std_devs;
		}

	public double[] getSS_dist_pca_mode_COV_min()
		{

			return dist_pca_mode_COV_min;
		}

	public double[] getSS_dist_pca_mode_COV_max()
		{

			return dist_pca_mode_COV_max;
		}

	public double[] getSS_dist_pca_mode_CORR_min()
		{

			return dist_pca_mode_CORR_min;
		}

	public double[] getSS_dist_pca_mode_CORR_max()
		{

			return dist_pca_mode_CORR_max;
		}

	public Matrix getDistance_matrix()
		{

			return distance_matrix;
		}

	public Matrix get_Subset_PDB_coordinates_dist()
		{

			return subset_PDB_coordinates_dist;
		}

	public Matrix getTransformed_subset_PDB_coordinates_dist()
		{

			return transformed_subset_PDB_coordinates_dist;
		}

	public Matrix getCov_dist()
		{

			return cov_dist;
		}

	public Matrix getCorr_dist()
		{

			return corr_dist;
		}

	public Matrix getTop_distance_evectors_COV()
		{

			return top_distance_evectors_COV;
		}

	public Matrix getTop_distance_evectors_CORR()
		{

			return top_distance_evectors_CORR;
		}

	public Matrix getNormed_projections_dist_COV()
		{

			return normed_projections_dist_COV;
		}

	public Matrix getProjections_dist_COV()
		{

			return projections_dist_COV;
		}

	public Matrix getProjections_dist_CORR()
		{

			return projections_dist_CORR;
		}

	public Matrix getNormed_projections_dist_CORR()
		{

			return normed_projections_dist_CORR;
		}

	public Matrix getWeighted_normed_projections_dist_COV()
		{

			return weighted_normed_projections_dist_COV;
		}

	public Matrix getWeighted_normed_projections_dist_CORR()
		{

			return weighted_normed_projections_dist_CORR;
		}

	public int getNumber_of_residues_dist()
		{

			return number_of_residues_dist;
		}

	public int getNumber_of_modes_dist()
		{

			return number_of_modes_dist;
		}

	public int[] getRes_list_dist()
		{

			return res_list_dist;
		}

	public int[] getRes_list_dist_orig()
		{

			return res_list_dist_orig;
		}

	public List<String> getChain_idents()
		{
			return chain_idents;
		}
}
