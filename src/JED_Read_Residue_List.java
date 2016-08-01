package jed;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import Jama.Matrix;

/**
 * JED class JED_Read_Residue_List: Reads a specified list of residues using the reference PDB file as a guide to select a subset.
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
 * @author Dr.Charles David
 */

public class JED_Read_Residue_List
{
	public String directory, description, pdb_ref_file, rl_SS;
	public int number_of_residues;
	public List<Integer> residue_list, residue_list_original, residues_read;
	public StringTokenizer sToken;
	public List<String> chain_idents, chain_idents_read, residue_ID_pairs_read;
	public int[] res_list, res_list_orig;
	public Matrix original_PDB_coordinates;

	/* ************************************************* CONSTRUCTORS **************************************************************************** */

	/**
	 * Constructor for reading residue lists
	 * 
	 * @param dir
	 *            The working directory
	 * @param desc
	 *            The job description
	 * @param res_list
	 *            The name of the residue list file
	 * @param pdb_ref
	 *            The PDB Reference File
	 * @param coords
	 *            The Coordinates Matrix (original_PDB_coordinates.txt)
	 */
	public JED_Read_Residue_List(String dir, String desc, String res_list, String pdb_ref, Matrix coords)
		{
			this.directory = dir;
			this.description = desc;
			this.rl_SS = res_list;
			this.pdb_ref_file = pdb_ref;
			this.original_PDB_coordinates = coords;
		}

	/* *************************************************** METHODS ******************************************************************************** */

	/**
	 * This method reads residue lists for Single Chain PDB files
	 */
	public void Read_Residue_List_Single()
		{
			try
				{
					JED_Read_PDBs rPDB = new JED_Read_PDBs(directory, description, pdb_ref_file);
					rPDB.read_Reference_PDB();
					residues_read = rPDB.get_Residues_Read();
					if (3 * residues_read.size() != original_PDB_coordinates.getRowDimension())
						{
							System.err.println("FATAL ERROR! The Reference PDB File DOES NOT MATCH the specified matrix of coordinates.");
							System.err.println("Terminating program execution.");
							System.exit(0);
						}

					File residues = new File(directory + rl_SS);
					BufferedReader residue_reader = new BufferedReader(new FileReader(residues));
					residue_list_original = new ArrayList<Integer>(); // preserves residue numbering for accessing the PDB file.
					residue_list = new ArrayList<Integer>(); // adjusts residue numbering to access X,Y,Z packing in the coordinates matrix.
					String line;
					while ((line = residue_reader.readLine()) != null)
						{
							sToken = new StringTokenizer(line);
							String r = sToken.nextToken();
							int res = Integer.parseInt(r);
							residue_list_original.add(res);
							int res_index = residues_read.indexOf(res);
							if (res_index == -1)
								{
									System.err.println("FATAL ERROR! Requested Residue DOES NOT EXIST in the Reference PDB File: " + res);
									System.err.println("Terminating program execution.");
									System.exit(0);
								}
							residue_list.add(res_index);
						}

					number_of_residues = residue_list.size();

					res_list_orig = new int[number_of_residues]; // INT array preserves residue numbering for accessing the PDB file.
					res_list = new int[number_of_residues]; // INT array adjusts residue numbering to access X,Y,Z packing in the coordinates matrix.

					for (int i = 0; i < number_of_residues; i++)
						{
							res_list[i] = residue_list.get(i);
							res_list_orig[i] = residue_list_original.get(i);
						}
					residue_reader.close();
				} catch (IOException io)
				{
					System.err.println("IOException thrown. Could not read the residue list file: " + directory + rl_SS);
					System.err.println("Terminating program execution.");
					io.printStackTrace();
					System.exit(0);
				}

		}

	/**
	 * This method reads residue lists for Multi Chain PDB files
	 */
	public void Read_Residue_List_Multi()
		{
			try
				{
					JED_Read_PDBs rPDB = new JED_Read_PDBs(directory, description, pdb_ref_file);
					rPDB.read_Reference_PDB();
					residues_read = rPDB.get_Residues_Read();
					chain_idents_read = rPDB.get_Chain_IDs_Read();
					residue_ID_pairs_read = rPDB.get_Residue_ID_Pairs();

					if (3 * residues_read.size() != original_PDB_coordinates.getRowDimension())
						{
							System.err.println("FATAL ERROR! The Reference PDB File DOES NOT MATCH the specified matrix of coordinates.");
							System.err.println("Terminating program execution.");
							System.exit(0);
						}

					File residues = new File(directory + rl_SS);
					BufferedReader residue_reader = new BufferedReader(new FileReader(residues));
					residue_list_original = new ArrayList<Integer>(); // preserves residue numbering for accessing the PDB file.
					residue_list = new ArrayList<Integer>(); // adjusts residue numbering to access the coordinates matrix.
					chain_idents = new ArrayList<String>(); // list of all chain identifiers specified
					String line, Chain_ID;
					while ((line = residue_reader.readLine()) != null)
						{
							sToken = new StringTokenizer(line);
							Chain_ID = sToken.nextToken();
							chain_idents.add(Chain_ID);
							int res = Integer.parseInt(sToken.nextToken());
							residue_list_original.add(res);
							String ID_Pair = Chain_ID + res;
							int res_index = residue_ID_pairs_read.indexOf(ID_Pair);
							if (res_index == -1)
								{
									System.err.println("FATAL ERROR! Requested Residue DOES NOT EXIST in the Reference PDB File: " + Chain_ID + res);
									System.err.println("Terminating program execution.");
									System.exit(0);
								}
							residue_list.add(res_index);
						}
					residue_reader.close();

					number_of_residues = residue_list.size();
					res_list_orig = new int[number_of_residues]; // INT array holing the original residue list values.
					res_list = new int[number_of_residues]; // INT array holing the residue list values adjusted for the matrix packing.

					for (int i = 0; i < number_of_residues; i++)
						{
							res_list[i] = residue_list.get(i);
							res_list_orig[i] = residue_list_original.get(i);
						}

				} catch (IOException io)
				{
					System.err.println("IOException thrown. Could not read the residue list file: " + directory + rl_SS);
					System.err.println("Terminating program execution.");
					io.printStackTrace();
					System.exit(0);
				}
		}

	/* *************************************************** GETTERS ******************************************************************************** */

	public int getNumber_of_residues()
		{
			return number_of_residues;
		}

	public int[] getRes_list()
		{
			return res_list;
		}

	public int[] getRes_list_orig()
		{
			return res_list_orig;
		}

	public List<Integer> getResidue_list()
		{
			return residue_list;
		}

	public List<Integer> getResidue_list_original()
		{
			return residue_list_original;
		}

	public List<Integer> getResidues_read()
		{
			return residues_read;
		}

	public List<String> getChain_idents()
		{
			return chain_idents;
		}

	public List<String> getChain_idents_read()
		{
			return chain_idents_read;
		}

	public List<String> getResidue_ID_pairs_read()
		{
			return residue_ID_pairs_read;
		}
}
