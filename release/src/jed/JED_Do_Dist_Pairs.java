package jed;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.List;
import java.util.StringTokenizer;

import Jama.Matrix;

/**
 * JED class JED_Do_Dist_Pairs: Top class for the Distance Pair analysis.
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
public class JED_Do_Dist_Pairs
{

	int number_of_Chains, number_of_pairs, reference_column, number_of_modes_dist_pairs;
	String directory, out_dir, description, pdb_ref_file, rl_dist_pairs, path, line;
	double trace_dist_COV, trace_dist_CORR, trace_dist_PCORR, cond_cov, cond_corr, cond_pcorr, det_cov, det_corr, det_pcorr, rank_cov, rank_corr, rank_pcorr, z_cutoff;
	List<Integer> residues_read, residue_list_dist1, residue_list_dist_orig1, residue_list_dist2, residue_list_dist_orig2;
	List<String> lines, chain_idents1, chain_idents2, chainIDs;
	List<Double> top_distance_eigenvalues_COV, top_distance_eigenvalues_CORR, top_distance_eigenvalues_PCORR;
	int[] res_list_dist1, res_list_dist_orig1, res_list_dist2, res_list_dist_orig2;
	double[] residue_distance_means, residue_distance_std_devs, dist_pca_mode_COV_min, dist_pca_mode_COV_max, dist_pca_mode_CORR_min, dist_pca_mode_CORR_max;
	Matrix distance_matrix, distance_matrix_cond, original_PDB_coordinates, subset_PDB_coordinates_dist, transformed_subset_PDB_coordinates_dist, cov_dist, corr_dist, top_distance_evectors_COV,
			top_distance_evectors_CORR, top_distance_evectors_PCORR, normed_projections_dist_COV, projections_dist_COV, projections_dist_CORR, normed_projections_dist_CORR,
			weighted_normed_projections_dist_COV, weighted_normed_projections_dist_CORR, Z_scores, counts;
	StringTokenizer sToken;
	NumberFormat nf;
	RoundingMode rm;

	/* ***************************************************** CONSTRUCTORS **************************************************************************** */

	/**
	 * Constructor that initiates the Distance Pair Analysis:
	 *
	 * @param dir
	 *            The working directory
	 * @param desc
	 *            The job description
	 * @param ref_file
	 *            The PDB Reference file
	 * @param res_dist_pairs
	 *            The list of Residue Pairs to analyze
	 * @param ref_col
	 *            The Reference Column in the Coordinates Matrix
	 * @param modes_dist_pairs
	 *            The number of Distance Pair PCA modes to process
	 * @param PDB_coordinates
	 *            The Matrix of Coordinates
	 */
	@SuppressWarnings("unused")
	public JED_Do_Dist_Pairs(String dir, String desc, String ref_file, String res_dist_pairs, int ref_col, int modes_dist_pairs, Matrix PDB_coordinates)
		{
			super();
			directory = dir;
			description = desc;
			pdb_ref_file = ref_file;
			out_dir = directory + "JED_RESULTS_" + description + "/dpPCA/";
			boolean exist = new File(out_dir).exists();
			if (!exist)
				{
					boolean success = (new File(out_dir)).mkdirs();
				}
			rl_dist_pairs = res_dist_pairs;
			reference_column = ref_col;
			number_of_modes_dist_pairs = modes_dist_pairs;
			original_PDB_coordinates = PDB_coordinates;

			nf = NumberFormat.getInstance();
			rm = RoundingMode.HALF_UP;
			nf.setRoundingMode(rm);
			nf.setMaximumFractionDigits(3);
			nf.setMinimumFractionDigits(3);
		}

	/* ************************************************** DRIVER METHODS ******************************************************************************** */

	/**
	 * Method for performing dPCA, processing Single Chain PDBs with no chain IDs.
	 */
	public void do_Dist()
		{
			read_Residues_Pairs_Single();
			get_Residue_Pair_Distances();
			remove_Outliers_Z_Score();
			get_Distance_Pair_PCA();
			get_Distance_Pair_DVPs();
			write_Distance_Pair_Stats();
		}

	/**
	 * Method for performing dPCA, processing Multi chain PDBs with chain IDs.
	 */
	public void do_Dist_Multi()
		{
			read_Residues_Pairs_Multi();
			get_Residue_Pair_Distances();
			remove_Outliers_Z_Score();
			get_Distance_Pair_PCA();
			get_Distance_Pair_DVPs();
			write_Distance_Pair_Stats_Multi();
		}

	/* ******************************************************** METHODS ******************************************************************************** */

	private void get_Distance_Pair_DVPs()
		{
			JED_Get_Distance_Pair_DVPs SS_dvp_COV = new JED_Get_Distance_Pair_DVPs(distance_matrix, top_distance_evectors_COV, top_distance_eigenvalues_COV, reference_column, directory, description,
					"COV", number_of_modes_dist_pairs);
				{
					// COVARIANCE METHOD
					SS_dvp_COV.get_Distance_Pair_DV_Series();
				}
			JED_Get_Distance_Pair_DVPs SS_dvp_CORR = new JED_Get_Distance_Pair_DVPs(distance_matrix, top_distance_evectors_CORR, top_distance_eigenvalues_CORR, reference_column, directory,
					description, "CORR", number_of_modes_dist_pairs);
				{
					// CORRELATION METHOD
					SS_dvp_CORR.get_Distance_Pair_DV_Series();
				}
			JED_Get_Distance_Pair_DVPs SS_dvp_PCORR = new JED_Get_Distance_Pair_DVPs(distance_matrix, top_distance_evectors_CORR, top_distance_eigenvalues_CORR, reference_column, directory,
					description, "PCORR", number_of_modes_dist_pairs);
				{
					// PARTIAL CORRELATION METHOD
					SS_dvp_PCORR.get_Distance_Pair_DV_Series();
				}
		}

	private void get_Distance_Pair_PCA()
		{
			JED_Get_Distance_Pair_PCA dp_pca = new JED_Get_Distance_Pair_PCA(distance_matrix_cond, directory, description, number_of_modes_dist_pairs, number_of_pairs);
				{
					dp_pca.get_Distance_Pair_PCA();

					residue_distance_means = dp_pca.getResidue_means().getColumnPackedCopy();
					residue_distance_std_devs = dp_pca.getResidue_std_devs().getColumnPackedCopy();

					// COVARIANCE METHOD
					trace_dist_COV = dp_pca.getTrace_COV();
					cond_cov = dp_pca.get_cond_COV();
					det_cov = dp_pca.getDet_COV();
					rank_cov = dp_pca.getRank_COV();
					top_distance_evectors_COV = dp_pca.getTop_evectors_dist_COV();
					top_distance_eigenvalues_COV = dp_pca.getTop_eigenvalues_COV();

					// CORRELATION METHOD
					trace_dist_CORR = dp_pca.getTrace_CORR();
					top_distance_evectors_CORR = dp_pca.getTop_evectors_dist_CORR();
					top_distance_eigenvalues_CORR = dp_pca.getTop_eigenvalues_CORR();

					// PARTIAL CORRELATION METHOD
					trace_dist_PCORR = dp_pca.getTrace_PCORR();
					top_distance_evectors_PCORR = dp_pca.getTop_evectors_PCORR();
					top_distance_eigenvalues_PCORR = dp_pca.getTop_eigenvalues_PCORR();
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
							Matrix_IO.write_Matrix(Z_scores, out_dir + "ss_" + number_of_pairs + "_Residue_Pairs_Distance_Z_scores.txt", 6, 1);
							counts = cv.get_counts();
							Matrix_IO.write_Matrix(counts, out_dir + "ss_" + number_of_pairs + "_Residue_Pairs_Outliers_Per_Variable.txt", 6, 0);
						}
				}
			if (z_cutoff == 0) distance_matrix_cond = distance_matrix;
		}

	private void write_Distance_Pair_Stats()
		{
			try
				{
					path = out_dir + "ss_" + number_of_pairs + "_Residue_Pairs_Distance_Stats.txt";
					File d_stats = new File(path);
					BufferedWriter d_stats_writer = new BufferedWriter(new FileWriter(d_stats));
					d_stats_writer.write("MEANs and STANDARD DEVIATIONS for the Residue Pair Distances: " + "\n");
					d_stats_writer.write(String.format("%-12s%-16s%-16s%-16s%n", "Res1", "Res2", "Mean", "Std_Dev"));
					for (int i = 0; i < number_of_pairs; i++)
						{
								{
									d_stats_writer.write(String.format("%-12s%-16s%-16s%-16s%n", res_list_dist_orig1[i], res_list_dist_orig2[i], nf.format(residue_distance_means[i]),
											nf.format(residue_distance_std_devs[i])));
								}
						}
					d_stats_writer.close();

				} catch (IOException io)
				{
					System.err.println("IOException thrown. Could not write the file: " + path);
					System.err.println("Terminating program execution.");
					io.printStackTrace();
					System.exit(0);
				}
		}

	private void get_Residue_Pair_Distances()
		{
			JED_Get_Distances_for_Residue_Pairs dist_pairs = new JED_Get_Distances_for_Residue_Pairs(original_PDB_coordinates, res_list_dist1, res_list_dist2, directory, description);
				{
					distance_matrix = dist_pairs.Get_Distance_Pairs();
				}
		}

	private void write_Distance_Pair_Stats_Multi()
		{
			try
				{
					path = out_dir + "ss_" + number_of_pairs + "_Residue_Pairs_Distance_Stats.txt";
					File d_stats = new File(path);
					BufferedWriter d_stats_writer = new BufferedWriter(new FileWriter(d_stats));
					d_stats_writer.write("MEANs and STANDARD DEVIATIONS for the Residue Pair Distances: " + "\n");
					d_stats_writer.write(String.format("%-12s%-16s%-16s%-16s%n", "Res1", "Res2", "Mean", "Std_Dev"));
					for (int i = 0; i < number_of_pairs; i++)
						{
								{
									d_stats_writer.write(String.format("%-12s%-16s%-16s%-16s%n", chain_idents1.get(i) + res_list_dist_orig1[i], chain_idents2.get(i) + res_list_dist_orig2[i],
											nf.format(residue_distance_means[i]), nf.format(residue_distance_std_devs[i])));
								}
						}
					d_stats_writer.close();

				} catch (IOException io)
				{
					System.out.println("IOException thrown. Could not write the file: " + path);
					System.err.println("Terminating program execution.");
					io.printStackTrace();
					System.exit(0);
				}
		}

	private void read_Residues_Pairs_Single()
		{
			JED_Read_Residue_Pair_List rpl = new JED_Read_Residue_Pair_List(directory, description, rl_dist_pairs, pdb_ref_file, original_PDB_coordinates, number_of_modes_dist_pairs);
			rpl.Read_Residue_List_Pairs_Single();

			residue_list_dist1 = rpl.getResidue_list1();
			residue_list_dist_orig1 = rpl.getResidue_list_original1();
			residue_list_dist2 = rpl.getResidue_list2();
			residue_list_dist_orig2 = rpl.getResidue_list_original2();
			res_list_dist1 = rpl.getRes_list1();
			res_list_dist_orig1 = rpl.getRes_list_orig1();
			res_list_dist2 = rpl.getRes_list2();
			res_list_dist_orig2 = rpl.getRes_list_orig2();
			number_of_pairs = residue_list_dist1.size();

			if (number_of_modes_dist_pairs > number_of_pairs)
				{
					System.err.println("ERROR!");
					System.err.println("Number of Distance Pair Modes REQUESTED: " + number_of_modes_dist_pairs);
					System.err.println("Number of Distance Pair Modes AVAILABLE: " + number_of_pairs);
					System.err.println("Possible Number of Distance Pair Modes is ALWAYS <= Number of Residues Pairs.");
					System.err.println("Setting Number of Distance Pair  Modes to Number of Residues Pairs:");
					number_of_modes_dist_pairs = number_of_pairs;
				}
		}

	private void read_Residues_Pairs_Multi()
		{
			JED_Read_Residue_Pair_List rpl = new JED_Read_Residue_Pair_List(directory, description, rl_dist_pairs, pdb_ref_file, original_PDB_coordinates, number_of_modes_dist_pairs);
			rpl.Read_Residue_List_Pairs_Multi();

			residue_list_dist1 = rpl.getResidue_list1();
			residue_list_dist_orig1 = rpl.getResidue_list_original1();
			residue_list_dist2 = rpl.getResidue_list2();
			residue_list_dist_orig2 = rpl.getResidue_list_original2();
			chain_idents1 = rpl.getChain_idents1();
			chain_idents2 = rpl.getChain_idents2();
			res_list_dist1 = rpl.getRes_list1();
			res_list_dist_orig1 = rpl.getRes_list_orig1();
			res_list_dist2 = rpl.getRes_list2();
			res_list_dist_orig2 = rpl.getRes_list_orig2();
			number_of_pairs = residue_list_dist1.size();

			if (number_of_modes_dist_pairs > number_of_pairs)
				{
					System.err.println("ERROR!");
					System.err.println("Number of Distance Pair Modes REQUESTED: " + number_of_modes_dist_pairs);
					System.err.println("Number of Distance Pair Modes AVAILABLE: " + number_of_pairs);
					System.err.println("Possible Number of Distance Pair Modes is ALWAYS <= Number of Residues Pairs.");
					System.err.println("Setting Number of Distance Pair  Modes to Number of Residues Pairs:");
					number_of_modes_dist_pairs = number_of_pairs;
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

	/**
	 * @return the number_of_Chains
	 */
	public int getNumber_of_Chains()
		{
			return number_of_Chains;
		}

	/**
	 * @return the number_of_pairs
	 */
	public int getNumber_of_pairs()
		{
			return number_of_pairs;
		}

	/**
	 * @return the reference_column
	 */
	public int getReference_column()
		{
			return reference_column;
		}

	/**
	 * @return the number_of_modes_dist_pairs
	 */
	public int getNumber_of_modes_dist_pairs()
		{
			return number_of_modes_dist_pairs;
		}

	/**
	 * @return the directory
	 */
	public String getDirectory()
		{
			return directory;
		}

	/**
	 * @return the out_dir
	 */
	public String getOut_dir()
		{
			return out_dir;
		}

	/**
	 * @return the description
	 */
	public String getDescription()
		{
			return description;
		}

	/**
	 * @return the pdb_ref_file
	 */
	public String getPdb_ref_file()
		{
			return pdb_ref_file;
		}

	/**
	 * @return the rl_dist_pairs
	 */
	public String getRl_dist_pairs()
		{
			return rl_dist_pairs;
		}

	/**
	 * @return the path
	 */
	public String getPath()
		{
			return path;
		}

	/**
	 * @return the line
	 */
	public String getLine()
		{
			return line;
		}

	/**
	 * @return the trace_dist_COV
	 */
	public double getTrace_dist_COV()
		{
			return trace_dist_COV;
		}

	/**
	 * @return the trace_dist_CORR
	 */
	public double getTrace_dist_CORR()
		{
			return trace_dist_CORR;
		}

	/**
	 * @return the trace_dist_PCORR
	 */
	public double getTrace_dist_PCORR()
		{
			return trace_dist_PCORR;
		}

	/**
	 * @return the cond_cov
	 */
	public double getCond_cov()
		{
			return cond_cov;
		}

	/**
	 * @return the cond_corr
	 */
	public double getCond_corr()
		{
			return cond_corr;
		}

	/**
	 * @return the cond_pcorr
	 */
	public double getCond_pcorr()
		{
			return cond_pcorr;
		}

	/**
	 * @return the det_cov
	 */
	public double getDet_cov()
		{
			return det_cov;
		}

	/**
	 * @return the det_corr
	 */
	public double getDet_corr()
		{
			return det_corr;
		}

	/**
	 * @return the det_pcorr
	 */
	public double getDet_pcorr()
		{
			return det_pcorr;
		}

	/**
	 * @return the rank_cov
	 */
	public double getRank_cov()
		{
			return rank_cov;
		}

	/**
	 * @return the rank_corr
	 */
	public double getRank_corr()
		{
			return rank_corr;
		}

	/**
	 * @return the rank_pcorr
	 */
	public double getRank_pcorr()
		{
			return rank_pcorr;
		}

	/**
	 * @return the z_cutoff
	 */
	public double getZ_cutoff()
		{
			return z_cutoff;
		}

	/**
	 * @return the residues_read
	 */
	public List<Integer> getResidues_read()
		{
			return residues_read;
		}

	/**
	 * @return the residue_list_dist1
	 */
	public List<Integer> getResidue_list_dist1()
		{
			return residue_list_dist1;
		}

	/**
	 * @return the residue_list_dist_orig1
	 */
	public List<Integer> getResidue_list_dist_orig1()
		{
			return residue_list_dist_orig1;
		}

	/**
	 * @return the residue_list_dist2
	 */
	public List<Integer> getResidue_list_dist2()
		{
			return residue_list_dist2;
		}

	/**
	 * @return the residue_list_dist_orig2
	 */
	public List<Integer> getResidue_list_dist_orig2()
		{
			return residue_list_dist_orig2;
		}

	/**
	 * @return the lines
	 */
	public List<String> getLines()
		{
			return lines;
		}

	/**
	 * @return the chain_idents1
	 */
	public List<String> getChain_idents1()
		{
			return chain_idents1;
		}

	/**
	 * @return the chain_idents2
	 */
	public List<String> getChain_idents2()
		{
			return chain_idents2;
		}

	/**
	 * @return the chainIDs
	 */
	public List<String> getChainIDs()
		{
			return chainIDs;
		}

	/**
	 * @return the top_distance_eigenvalues_COV
	 */
	public List<Double> getTop_distance_eigenvalues_COV()
		{
			return top_distance_eigenvalues_COV;
		}

	/**
	 * @return the top_distance_eigenvalues_CORR
	 */
	public List<Double> getTop_distance_eigenvalues_CORR()
		{
			return top_distance_eigenvalues_CORR;
		}

	/**
	 * @return the top_distance_eigenvalues_PCORR
	 */
	public List<Double> getTop_distance_eigenvalues_PCORR()
		{
			return top_distance_eigenvalues_PCORR;
		}

	/**
	 * @return the res_list_dist1
	 */
	public int[] getRes_list_dist1()
		{
			return res_list_dist1;
		}

	/**
	 * @return the res_list_dist_orig1
	 */
	public int[] getRes_list_dist_orig1()
		{
			return res_list_dist_orig1;
		}

	/**
	 * @return the res_list_dist2
	 */
	public int[] getRes_list_dist2()
		{
			return res_list_dist2;
		}

	/**
	 * @return the res_list_dist_orig2
	 */
	public int[] getRes_list_dist_orig2()
		{
			return res_list_dist_orig2;
		}

	/**
	 * @return the residue_distance_means
	 */
	public double[] getResidue_distance_means()
		{
			return residue_distance_means;
		}

	/**
	 * @return the residue_distance_std_devs
	 */
	public double[] getResidue_distance_std_devs()
		{
			return residue_distance_std_devs;
		}

	/**
	 * @return the dist_pca_mode_COV_min
	 */
	public double[] getDist_pca_mode_COV_min()
		{
			return dist_pca_mode_COV_min;
		}

	/**
	 * @return the dist_pca_mode_COV_max
	 */
	public double[] getDist_pca_mode_COV_max()
		{
			return dist_pca_mode_COV_max;
		}

	/**
	 * @return the dist_pca_mode_CORR_min
	 */
	public double[] getDist_pca_mode_CORR_min()
		{
			return dist_pca_mode_CORR_min;
		}

	/**
	 * @return the dist_pca_mode_CORR_max
	 */
	public double[] getDist_pca_mode_CORR_max()
		{
			return dist_pca_mode_CORR_max;
		}

	/**
	 * @return the distance_matrix
	 */
	public Matrix getDistance_matrix()
		{
			return distance_matrix;
		}

	/**
	 * @return the distance_matrix_cond
	 */
	public Matrix getDistance_matrix_cond()
		{
			return distance_matrix_cond;
		}

	/**
	 * @return the original_PDB_coordinates
	 */
	public Matrix getOriginal_PDB_coordinates()
		{
			return original_PDB_coordinates;
		}

	/**
	 * @return the subset_PDB_coordinates_dist
	 */
	public Matrix getSubset_PDB_coordinates_dist()
		{
			return subset_PDB_coordinates_dist;
		}

	/**
	 * @return the transformed_subset_PDB_coordinates_dist
	 */
	public Matrix getTransformed_subset_PDB_coordinates_dist()
		{
			return transformed_subset_PDB_coordinates_dist;
		}

	/**
	 * @return the cov_dist
	 */
	public Matrix getCov_dist()
		{
			return cov_dist;
		}

	/**
	 * @return the corr_dist
	 */
	public Matrix getCorr_dist()
		{
			return corr_dist;
		}

	/**
	 * @return the top_distance_evectors_COV
	 */
	public Matrix getTop_distance_evectors_COV()
		{
			return top_distance_evectors_COV;
		}

	/**
	 * @return the top_distance_evectors_CORR
	 */
	public Matrix getTop_distance_evectors_CORR()
		{
			return top_distance_evectors_CORR;
		}

	/**
	 * @return the top_distance_evectors_PCORR
	 */
	public Matrix getTop_distance_evectors_PCORR()
		{
			return top_distance_evectors_PCORR;
		}

	/**
	 * @return the normed_projections_dist_COV
	 */
	public Matrix getNormed_projections_dist_COV()
		{
			return normed_projections_dist_COV;
		}

	/**
	 * @return the projections_dist_COV
	 */
	public Matrix getProjections_dist_COV()
		{
			return projections_dist_COV;
		}

	/**
	 * @return the projections_dist_CORR
	 */
	public Matrix getProjections_dist_CORR()
		{
			return projections_dist_CORR;
		}

	/**
	 * @return the normed_projections_dist_CORR
	 */
	public Matrix getNormed_projections_dist_CORR()
		{
			return normed_projections_dist_CORR;
		}

	/**
	 * @return the weighted_normed_projections_dist_COV
	 */
	public Matrix getWeighted_normed_projections_dist_COV()
		{
			return weighted_normed_projections_dist_COV;
		}

	/**
	 * @return the weighted_normed_projections_dist_CORR
	 */
	public Matrix getWeighted_normed_projections_dist_CORR()
		{
			return weighted_normed_projections_dist_CORR;
		}

	/**
	 * @return the z_scores
	 */
	public Matrix getZ_scores()
		{
			return Z_scores;
		}

	/**
	 * @return the counts
	 */
	public Matrix getCounts()
		{
			return counts;
		}

	/**
	 * @return the number format
	 */
	public NumberFormat getNf()
		{
			return nf;
		}

	/**
	 * @return the rounding
	 */
	public RoundingMode getRm()
		{
			return rm;
		}
}
