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
 * JED class JED_Read_Residue_Pair_List: Reads a specified list of residue pairs using the reference PDB file as a guide to select a subset.
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

public class JED_Read_Residue_Pair_List
	{
		public String directory, description, pdb_ref_file, rl_pairs;
		public int number_of_pairs, number_of_modes;
		public List<Integer> residue_list1, residue_list2, residue_list_original1, residue_list_original2, residues_read;
		public StringTokenizer sToken;
		public List<String> chain_idents1, chain_idents2, chain_idents_read, residue_ID_pairs_read;
		public int[] res_list1, res_list2, res_list_orig1, res_list_orig2;
		public Matrix original_PDB_coordinates;

		/* ************************************************* CONSTRUCTORS **************************************************************************** */

		/**
		 * Constructor to read a residue pair list
		 *
		 * @param dir
		 *            The working directory
		 * @param desc
		 *            The job description
		 * @param pair_list
		 *            The name of the residue pair list
		 * @param pdb_ref
		 *            The PDB Reference File
		 * @param coords
		 *            The coordinates matrix (original_PDB_coordinates.txt)
		 * @param num_modes
		 *            The number of dpPCA modes to process
		 */
		public JED_Read_Residue_Pair_List(String dir, String desc, String pair_list, String pdb_ref, Matrix coords, int num_modes)
			{
				this.directory = dir;
				this.description = desc;
				this.rl_pairs = pair_list;
				this.pdb_ref_file = pdb_ref;
				this.original_PDB_coordinates = coords;
				this.number_of_modes = num_modes;
			}

		/* *************************************************** METHODS ******************************************************************************** */

		/**
		 * This method reads the residue pair list for Single Chain PDB files.
		 */
		public void Read_Residue_List_Pairs_Single()
			{
				try
					{
						JED_Read_PDBs rPDB = new JED_Read_PDBs(directory, description, pdb_ref_file);
						rPDB.read_Reference_PDB();
						residues_read = rPDB.get_Residues_Read();
						if (3 * residues_read.size() != original_PDB_coordinates.getRowDimension())
							{
								System.err.println("ERROR! The Reference PDB File DOES NOT MATCH the specified matrix of coordinates.");
								System.err.println("Terminating program execution.");
								System.exit(0);
							}

						File residues = new File(directory + rl_pairs);
						BufferedReader residue_reader = new BufferedReader(new FileReader(residues));

						residue_list1 = new ArrayList<Integer>();
						residue_list_original1 = new ArrayList<Integer>();
						residue_list2 = new ArrayList<Integer>();
						residue_list_original2 = new ArrayList<Integer>();

						String line, element;

						while ((line = residue_reader.readLine()) != null)
							{
								sToken = new StringTokenizer(line);
								if (sToken.countTokens() < 2)
									{
										System.err.println("ERROR! Single Chain PDB residue pair list must have 2 columns!");
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								element = sToken.nextToken();
								int res1 = Integer.parseInt(element);
								residue_list_original1.add(res1);
								int res1_index = residues_read.indexOf(res1);
								residue_list1.add(res1_index);
								element = sToken.nextToken();
								if (sToken.hasMoreTokens())
									System.err.println("ERROR! Residue list for SINGLE chain PDB file should have a ONLY 2 columns of numbers!");
								int res2 = Integer.parseInt(element);
								residue_list_original2.add(res2);
								int res2_index = residues_read.indexOf(res2);
								residue_list2.add(res2_index);
								if (res1_index == -1 || res2_index == -1)
									{
										System.err.println("ERROR! Requested Residue Pair DOES NOT EXIST in the Reference PDB File: " + res1 + "\t" + res2);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
							}
						residue_reader.close();
						number_of_pairs = residue_list1.size();

						if (number_of_modes > (number_of_pairs))
							{
								System.err.println("ERROR!");
								System.err.println("Number of Distance Pair Modes REQUESTED: " + number_of_modes);
								System.err.println("Number of Distance Pair Modes AVAILABLE: " + number_of_pairs);
								System.err.println("Possible number of Distance Pair Modes is ALWAYS <= Number of Residues Pairs.");
								System.err.println("Terminating program execution.");
								System.exit(0);
							}

						res_list1 = new int[number_of_pairs];
						res_list_orig1 = new int[number_of_pairs];
						res_list2 = new int[number_of_pairs];
						res_list_orig2 = new int[number_of_pairs];

						for (int i = 0; i < number_of_pairs; i++)
							{
								res_list1[i] = residue_list1.get(i);
								res_list_orig1[i] = residue_list_original1.get(i);
								res_list2[i] = residue_list2.get(i);
								res_list_orig2[i] = residue_list_original2.get(i);
							}
					} catch (IOException io)
					{
						System.out.println("IOException thrown. Could not read the residue pair file: " + directory + rl_pairs);
						System.err.println("Terminating program execution.");
						io.printStackTrace();
						System.exit(0);
					}

			}

		/**
		 * This method reads the residue pair list for Multi Chain PDB files.
		 */
		public void Read_Residue_List_Pairs_Multi()
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
								System.err.println("ERROR! The Reference PDB File DOES NOT MATCH the specified matrix of coordinates.");
								System.err.println("Terminating program execution.");
								System.exit(0);
							}

						File residues = new File(directory + rl_pairs);
						BufferedReader residue_reader = new BufferedReader(new FileReader(residues));

						residue_list1 = new ArrayList<Integer>();
						residue_list_original1 = new ArrayList<Integer>();
						residue_list2 = new ArrayList<Integer>();
						residue_list_original2 = new ArrayList<Integer>();
						chain_idents1 = new ArrayList<String>();
						chain_idents2 = new ArrayList<String>();

						String line = null, element;
						while ((line = residue_reader.readLine()) != null)
							{
								// Reading the first list of residues: Column 1: chain IDs and Column 2: residue numbers

								sToken = new StringTokenizer(line);
								if (sToken.countTokens() < 4)
									{
										System.err.println("ERROR! Multi Chain PDB residue pairs list must have 4 columns!");
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								String chain_ID1 = sToken.nextToken();
								chain_idents1.add(chain_ID1);
								element = sToken.nextToken();
								int res1 = Integer.parseInt(element);
								residue_list_original1.add(res1);
								String ID_Pair1 = chain_ID1 + res1;
								int res_index1 = residue_ID_pairs_read.indexOf(ID_Pair1);
								if (res_index1 == -1)
									{
										System.err.println("ERROR! Requested Residue DOES NOT EXIST in the Reference PDB File: " + chain_ID1 + res1);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								residue_list1.add(res_index1);

								// Reading the second list of residues: Column 3: chain IDs and Column 4: residue numbers

								String chain_ID2 = sToken.nextToken();
								chain_idents2.add(chain_ID2);
								element = sToken.nextToken();
								if (sToken.hasMoreTokens())
									System.err.println("ERROR! MULTI chain PDB residue pair list should only have 4 columns: 2 Res# and 2 ChainID!");
								int res2 = Integer.parseInt(element);
								residue_list_original2.add(res2);
								String ID_Pair2 = chain_ID2 + res2;
								int res_index2 = residue_ID_pairs_read.indexOf(ID_Pair2);
								if (res_index2 == -1)
									{
										System.err.println("ERROR! Requested Residue DOES NOT EXIST in the Reference PDB File: " + chain_ID2 + res2);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								residue_list2.add(res_index2);
							}
						residue_reader.close();
						number_of_pairs = residue_list1.size();

						if (number_of_modes > (number_of_pairs))
							{
								System.err.println("ERROR!");
								System.err.println("Number of Distance Pair Modes REQUESTED: " + number_of_modes);
								System.err.println("Number of Distance Pair Modes AVAILABLE: " + number_of_pairs);
								System.err.println("Possible number of Distance Pair Modes is AWAYS <= Number of Residue Pairs.");
								System.err.println("Terminating program execution.");
								System.exit(0);
							}

						res_list1 = new int[number_of_pairs];
						res_list_orig1 = new int[number_of_pairs];
						res_list2 = new int[number_of_pairs];
						res_list_orig2 = new int[number_of_pairs];

						for (int i = 0; i < number_of_pairs; i++)
							{
								res_list1[i] = residue_list1.get(i);
								res_list_orig1[i] = residue_list_original1.get(i);
								res_list2[i] = residue_list2.get(i);
								res_list_orig2[i] = residue_list_original2.get(i);
							}

					} catch (IOException io)
					{
						System.out.println("IOException thrown. Could not read the residue list file: " + directory + rl_pairs);
						io.printStackTrace();
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
			}

		/* *************************************************** GETTERS ******************************************************************************** */

		public List<Integer> getResidue_list1()
			{
				return residue_list1;
			}

		public List<Integer> getResidue_list2()
			{
				return residue_list2;
			}

		public List<Integer> getResidue_list_original1()
			{
				return residue_list_original1;
			}

		public List<Integer> getResidue_list_original2()
			{
				return residue_list_original2;
			}

		public List<Integer> getResidues_read()
			{
				return residues_read;
			}

		public List<String> getChain_idents1()
			{
				return chain_idents1;
			}

		public List<String> getChain_idents2()
			{
				return chain_idents2;
			}

		public List<String> getChain_idents_read()
			{
				return chain_idents_read;
			}

		public List<String> getResidue_ID_pairs_read()
			{
				return residue_ID_pairs_read;
			}

		public int[] getRes_list1()
			{
				return res_list1;
			}

		public int[] getRes_list2()
			{
				return res_list2;
			}

		public int[] getRes_list_orig1()
			{
				return res_list_orig1;
			}

		public int[] getRes_list_orig2()
			{
				return res_list_orig2;
			}

	}
