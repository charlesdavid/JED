package jed;

import java.io.File;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

import Jama.Matrix;

/**
 * JED class JED_Do_NO_PCA. Top class for the Pre-Processing run the creates the Coordinates Matrix.
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
public class JED_Do_No_PCA
{

	public int number_of_residues, reference_column;
	public String directory, description, out_dir, pdb_ref_file, rl;
	public double z_cutoff, percent;
	public List<Integer> all_residue_list, residue_count_offsets, chain_lengths;
	public List<String> chainID_list;
	public List<Double> conformation_rmsds, residue_rmsd_list, Z_Scores;
	public Matrix original_PDB_coordinates, transformed_subset_PDB_coordinates, trimmed_PDB_coords_COLS, adjusted_PDB_coordinates_ROWS, conf_Z_scores;
	public Vector<Atom> atoms;
	public StringTokenizer sToken;

	/* ************************************** CONSTRUCTORS ******************************************************************************** */

	/**
	 * Initiates the Pre-Processing (NO PCA) Analysis for Single Chain PDB files:
	 * 
	 * @param dir
	 *            The working directory
	 * @param rl_name
	 *            The name of the residue list file
	 * @param pdb_res_read
	 *            The list of all residue numbers found in the PDB files
	 * @param desc
	 *            The job description
	 * @param ref_pdb
	 *            The PDB Reference file
	 * @param coordinates
	 *            The Matrix of Coordinates
	 * @param ref_col
	 *            The Reference Column in the Coordinates Matrix
	 */
	public JED_Do_No_PCA(String dir, String rl_name, List<Integer> pdb_res_read, String desc, String ref_pdb, Matrix coordinates, int ref_col)
		{

			super();
			this.directory = dir;
			this.rl = rl_name;
			this.all_residue_list = pdb_res_read;
			this.description = desc;
			this.pdb_ref_file = ref_pdb;
			this.original_PDB_coordinates = coordinates;
			this.reference_column = ref_col;
			this.out_dir = directory + "JED_RESULTS_" + description + "/";
			boolean exist = new File(out_dir).exists();
			if (!exist)
				{
					boolean success = (new File(out_dir)).mkdirs();
					if (!success)
						{
							System.err.println("Failed to create the output directory. Terminating execution.");
							System.exit(0);
						}
				}
		}

	/**
	 * Initiates the Pre-Processing (NO PCA) Analysis for Multi Chain PDB files:
	 * 
	 * @param dir
	 *            The working directory
	 * @param rl_name
	 *            The name of the residue list file
	 * @param chain_ids_read
	 *            List of all chain IDs found in the PDB files
	 * @param pdb_res_read
	 *            The list of all residue numbers found in the PDB files
	 * @param desc
	 *            The job description
	 * @param ref_pdb
	 *            The PDB Reference file
	 * @param coordinates
	 *            The Matrix of Coordinates
	 * @param ref_col
	 *            The Reference Column in the Coordinates Matrix
	 */
	public JED_Do_No_PCA(String dir, String rl_name, List<String> chain_ids_read, List<Integer> pdb_res_read, String desc, String ref_pdb, Matrix coordinates, int ref_col)
		{

			super();
			this.directory = dir;
			this.rl = rl_name;
			this.chainID_list = chain_ids_read;
			this.all_residue_list = pdb_res_read;
			this.description = desc;
			this.pdb_ref_file = ref_pdb;
			this.original_PDB_coordinates = coordinates;
			this.reference_column = ref_col;
			this.out_dir = directory + "JED_RESULTS_" + description + "/";
			boolean exist = new File(out_dir).exists();
			if (!exist)
				{
					boolean success = (new File(out_dir)).mkdirs();
					if (!success)
						{
							System.err.println("Failed to create the output directory. Terminating execution.");
							System.exit(0);
						}
				}
		}

	// ***************************************** SETTERS ***************************************************//

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

	/**
	 * Sets the Percent for outlier processing.
	 * 
	 * @param p
	 *            The Z cutoff
	 */
	public void set_percent_cutoff(double p)
		{

			percent = p;
		}

	// ***************************************** METHODS ***************************************************//

	/**
	 * Runs the NO PCA Methods for Single Chain PDB files
	 */
	public void do_No_PCA()
		{
			get_Transformed_Coords();
			JED_Edited_Subset_PDB pdb = new JED_Edited_Subset_PDB(directory, pdb_ref_file, residue_rmsd_list, all_residue_list);
				{
					pdb.set_Output_Directory(out_dir);
					pdb.do_edit_pdb();
					atoms = pdb.getAtoms();
				}
		}

	private void get_Transformed_Coords()
		{
			JED_Get_Transformed_Coordinates tf_coords = new JED_Get_Transformed_Coordinates(original_PDB_coordinates, reference_column, directory, description);
				{
					tf_coords.set_percent_cutoff(percent);
					tf_coords.set_z_cutoff(z_cutoff);

					transformed_subset_PDB_coordinates = tf_coords.get_SS_Transformed_coords();
					trimmed_PDB_coords_COLS = tf_coords.get_SS_transformed_coordinates_trimmed_COLS();
					adjusted_PDB_coordinates_ROWS = tf_coords.get_SS_transformed_coordinates_adjusted_ROWS();
					conformation_rmsds = tf_coords.get_SS_Conformation_RMSDs();
					residue_rmsd_list = tf_coords.get_SS_Residue_RMSDs();
					conf_Z_scores = tf_coords.get_conf_Z_scores();
					number_of_residues = tf_coords.getNumber_of_residues();
				}
		}

	/**
	 * Runs the NO PCA Methods for Multi Chain PDB files
	 */
	public void do_No_PCA_Multi()
		{
			get_Transformed_Coords();
			JED_Edited_Subset_PDB pdb = new JED_Edited_Subset_PDB(directory, pdb_ref_file, residue_rmsd_list, chainID_list, all_residue_list);
				{
					pdb.set_Output_Directory(out_dir);
					pdb.do_edit_pdb_multi();
					atoms = pdb.getAtoms();
				}
		}

	// ***************************************** GETTERS ***************************************************//

	public int getNumber_of_residues()
		{

			return number_of_residues;
		}

	public List<Integer> getResidue_list()
		{

			return all_residue_list;
		}

	public List<Double> get_conformation_rmsds()
		{

			return conformation_rmsds;
		}

	public List<Double> get_residue_rmsd_list()
		{

			return residue_rmsd_list;
		}

	public Matrix getTransformed_subset_PDB_coordinates()
		{

			return transformed_subset_PDB_coordinates;
		}

	public Vector<Atom> getAtoms()
		{

			return atoms;
		}

}
