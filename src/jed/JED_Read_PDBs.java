package jed;

import java.io.BufferedReader;
import java.io.File;
import java.util.List;
import java.util.Vector;

import Jama.Matrix;

/**
 * JED class JED_Read_PDBs: Top class for reading PDB files in a directory and getting the Coordinates Matrix.
 * This class also reads the PDB Reference File.
 * 
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
 * along with this program. If not, see <http://www.gnu.org/license>.
 * 
 * @author Dr. Charles David
 */

public class JED_Read_PDBs
	{

		String directory, out_dir, pdb_ref_file, description;
		Matrix original_reference_PDB_coordinates, original_PDB_coordinates;
		List<Integer> residues_read;
		List<String> pdb_file_names, chain_ids_read, residue_ID_pairs_read;
		Vector<Atom> atoms, atoms_CA;
		BufferedReader reader;
		FortranFormat formatter;
		PDB_File_Parser fp;
		JED_Get_Coordinates_from_PDBs coords;
		PDB_Coordinates data;
		boolean exist, success;

		// ***************************************** CONSTRUCTORS ***************************************************//

		/**
		 * Constructor that initiates the reading of PDB files
		 * 
		 * @param dir
		 *            The working directory
		 * @param desc
		 *            The job description
		 * @param ref_file
		 *            The PDB reference file
		 */
		public JED_Read_PDBs(String dir, String desc, String ref_file)
			{

				super();
				this.directory = dir;
				this.description = desc;
				out_dir = directory + "JED_RESULTS_" + description + "/";
				exist = new File(out_dir).exists();
				if (!exist) success = (new File(out_dir)).mkdirs();
				this.pdb_ref_file = ref_file;
			}

		// ***************************************** METHODS ***************************************************//

		/**
		 * This method reads all the PDB files in the working directory.
		 */
		public void read_PDBs()
			{
				coords = new JED_Get_Coordinates_from_PDBs(directory, description);
				original_PDB_coordinates = coords.get_Orig_PDB_Coords();
				pdb_file_names = coords.get_PDB_names();
				List_IO.write_String_List(pdb_file_names, out_dir + "PDB_READ_LOG.txt");
			}

		/**
		 * This method reads the PDB Reference File.
		 * The lists of all chain IDs and residue numbers found in the file are generated.
		 * Every atom in the file is stored in the vector of atoms.
		 */
		public void read_Reference_PDB()
			{
				PDB_IO pio = new PDB_IO(directory, pdb_ref_file);
				atoms = pio.Read_PDB();
				atoms_CA = pio.get_Atoms_CA();
				residues_read = pio.getResidues_read();
				chain_ids_read = pio.getChain_ids_read();
				residue_ID_pairs_read = pio.getResidue_id_pairs_read();
				original_reference_PDB_coordinates = PDB_Coordinates.get_PDB_Coords(atoms_CA);
			}

		// ***************************************** GETTERS ***************************************************//

		Matrix get_Original_PDB_Coordinates()
			{

				return original_PDB_coordinates;
			}

		/**
		 * @return the original_reference_PDB_coordinates
		 */
		public Matrix get_Original_Reference_PDB_coordinates()
			{
				return original_reference_PDB_coordinates;
			}

		List<String> get_PDB_File_Names()
			{

				return pdb_file_names;
			}

		public List<Integer> get_Residues_Read()
			{

				return residues_read;
			}

		public List<String> get_Chain_IDs_Read()
			{

				return chain_ids_read;
			}

		public List<String> get_Residue_ID_Pairs()
			{
				return residue_ID_pairs_read;
			}
	}
