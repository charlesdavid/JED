package jed;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Vector;

/**
 * JED class JED_Edited_Subset_PDB: Writes a subset of residues to a PDB file and replaces B-Factors with residue RMSF.
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
 * 
 * @author Dr. Charles David
 * 
 */
public class JED_Edited_Subset_PDB
{

	String directory, out_dir, pdb_data_file, CA = "CA";
	int number_of_residues;
	List<Double> residue_rmsd_list;
	List<Integer> residue_list;
	List<String> chainID_list;
	Vector<Atom> atoms;
	final double FLOOR = 1.00E-3, delta_y = 99;

	/* ************************************** CONSTRUCTORS ******************************************************************************** */

	/**
	 * Constructor for getting a subset of residues from a Single Chain PDB file.
	 * 
	 * @param res_list
	 *            The list of residues for the subset
	 * @param dir
	 *            The working directory
	 * @param ref_file
	 *            The PDB Reference File
	 */
	JED_Edited_Subset_PDB(List<Integer> res_list, String dir, String ref_file)

		{
			directory = dir;
			pdb_data_file = ref_file;
			residue_list = res_list;
			atoms = new Vector<Atom>();
			number_of_residues = residue_list.size();
		}

	/**
	 * Constructor for getting a subset of residues from a Multi Chain PDB file.
	 * 
	 * @param chain_ids
	 *            The list of Chain IDs of the residues for the subset
	 * @param res_list
	 *            The list of residues for the subset
	 * @param dir
	 *            The working directory
	 * @param ref_file
	 *            The PDB Reference File
	 */

	JED_Edited_Subset_PDB(List<String> chain_ids, List<Integer> res_list, String dir, String ref_file)

		{
			directory = dir;
			pdb_data_file = ref_file;
			residue_list = res_list;
			chainID_list = chain_ids;
			atoms = new Vector<Atom>();
			number_of_residues = residue_list.size();
		}

	/**
	 * Constructor for getting editing a subset of residues from a Single Chain PDB file:
	 * The B-Factors are replaced with the residue RMSFs
	 * 
	 * @param dir
	 *            The working directory
	 * @param ref_file
	 *            The PDB Reference File
	 * @param rmsds
	 *            The set of residue RMSDs (RMSFs)
	 * @param res
	 *            The list of residue numbers for the subset
	 */
	JED_Edited_Subset_PDB(String dir, String ref_file, List<Double> rmsds, List<Integer> res)

		{
			directory = dir;
			residue_rmsd_list = rmsds;
			pdb_data_file = ref_file;
			residue_list = res;
			atoms = new Vector<Atom>();
			number_of_residues = residue_rmsd_list.size();
		}

	/**
	 * Constructor for getting editing a subset of residues from a Single Chain PDB file:
	 * The B-Factors are replaced with the residue RMSFs
	 * 
	 * @param dir
	 *            The working directory
	 * @param ref_file
	 *            The PDB Reference File
	 * @param rmsds
	 *            The set of residue RMSDs (RMSFs)
	 * @param chain_ids
	 *            The list of chain IDs for the subset
	 * @param res
	 *            The list of residue numbers for the subset
	 */
	JED_Edited_Subset_PDB(String dir, String ref_file, List<Double> rmsds, List<String> chain_ids, List<Integer> res)
		{

			directory = dir;
			residue_rmsd_list = rmsds;
			pdb_data_file = ref_file;
			chainID_list = chain_ids;
			residue_list = res;
			atoms = new Vector<Atom>();
			number_of_residues = residue_rmsd_list.size();
		}

	/* ************************************** METHODS ******************************************************************************** */

	/**
	 * Writes the subset to a Single Chain PDB file
	 */
	public void get_pdb_subset()
		{

			JED_PDB_IO pdb_io = new JED_PDB_IO(directory, pdb_data_file);
			atoms = pdb_io.Read_PDB_Subset(residue_list);
			String ref_PDB_file_name = "ss_" + number_of_residues + ".pdb";
			pdb_io.Write_PDB(out_dir, ref_PDB_file_name, atoms);
		}

	/**
	 * Writes the subset to a Multi Chain PDB file
	 */
	public void get_pdb_subset_multi()
		{

			JED_PDB_IO pdb_io = new JED_PDB_IO(directory, pdb_data_file);
			atoms = pdb_io.Read_PDB_Subset_Multi(chainID_list, residue_list);
			String ref_PDB_file_name = "ss_" + number_of_residues + ".pdb";
			pdb_io.Write_PDB(out_dir, ref_PDB_file_name, atoms);
		}

	/**
	 * Writes the edited subset to a Single Chain PDB file
	 */
	public void do_edit_pdb()
		{

			JED_PDB_IO pdb_io = new JED_PDB_IO(directory, pdb_data_file);
			atoms = pdb_io.Read_PDB_Subset(residue_list);
			String ref_PDB_file_name = "ss_" + number_of_residues + "_RMSF_edited.pdb";
			edit_B_Factors_with_RMSDs();
			pdb_io.Write_PDB(out_dir, ref_PDB_file_name, atoms);
		}

	/**
	 * Writes the edited subset to a Multi Chain PDB file
	 */
	public void do_edit_pdb_multi()
		{

			JED_PDB_IO pdb_io = new JED_PDB_IO(directory, pdb_data_file);
			atoms = pdb_io.Read_PDB_Subset_Multi(chainID_list, residue_list);
			String ref_PDB_file_name = "ss_" + number_of_residues + "_RMSF_edited.pdb";
			edit_B_Factors_with_RMSDs();
			pdb_io.Write_PDB(out_dir, ref_PDB_file_name, atoms);
		}

	/**
	 * Replaces the B-Factors with the Residue RMSDs
	 */
	private void edit_B_Factors_with_RMSDs()
		{

			List<Double> sorted_res_rmsds = new ArrayList<Double>();
			sorted_res_rmsds.addAll(residue_rmsd_list);
			Collections.sort(sorted_res_rmsds, Collections.reverseOrder());
			double max_rmsd = sorted_res_rmsds.get(0);
			double min_rmsd = sorted_res_rmsds.get(residue_rmsd_list.size() - 1);
			double log_RR_max = Math.log10(max_rmsd);
			double log_RR_min = Math.log10(min_rmsd);
			double delta_x = ((log_RR_max) - (log_RR_min));
			double slope = (delta_y / delta_x);
			double y_min = (slope * log_RR_min);
			int d = 0;
			for (Atom a : atoms)
				{
					if (a.symbol.equals(CA))
						{
							double bff = (residue_rmsd_list.get(d));
							double log_bff = Math.log10(bff);
							double bf = ((slope * log_bff) - y_min);
							a.b_factor = bf;
							d++;
						} else
						{
							a.b_factor = 0;
						}
				}
		}

	/* ************************************** SETTERS ******************************************************************************** */

	/**
	 * Sets the output directory for this class.
	 * 
	 * @param dir
	 *            The output directory
	 */
	public void set_Output_Directory(String dir)
		{
			this.out_dir = dir;
		}

	/* ************************************** GETTERS ******************************************************************************** */

	/**
	 * Returns the Vector of atoms that comprise the subset.
	 * 
	 * @return
	 */
	public Vector<Atom> getAtoms()
		{

			return atoms;
		}
}
