package jed;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Vector;

/**
 * JED class PDB_IO: Handles the reading and writing of PDB files using a Fortran Formatter and a PDB File Parser.
 * Catches Exceptions thrown by the PDB File Parser class.
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

public class PDB_IO
{

	String read_directory, write_directory, read_name, write_name, path;
	Vector<Atom> atoms, atoms_CA, edited_atoms;
	List<String> chain_ids_read, residue_id_pairs_read;
	List<Integer> residues_read;
	File pdb, edited_pdb;
	BufferedReader pdb_reader;
	BufferedWriter pdb_writer;
	FortranFormat formatter;
	PDB_File_Parser parser;

	/* ************************** CONSTRUCTORS ******************************************************************* */

	/**
	 * Constructor for the PDB IO object:
	 * Instantiates the PDB reader and the Fortran Formatter.
	 * 
	 * @param dir
	 *            The read directory
	 * @param name
	 *            The read name
	 */
	public PDB_IO(String dir, String name)
		{

			read_name = name;
			read_directory = dir;
			pdb = new File(read_directory + read_name);
			try
				{
					pdb_reader = new BufferedReader(new FileReader(pdb));
				} catch (FileNotFoundException e)
				{
					System.err.println("Could not find the file: " + read_directory + read_name);
					e.printStackTrace();
				}
			formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
			formatter.setAddReturn(true);
			parser = new PDB_File_Parser();

		}
	
	public PDB_IO(String path)
	{

		this.path = path;
		pdb = new File(read_directory + read_name);
		try
			{
				pdb_reader = new BufferedReader(new FileReader(pdb));
			} catch (FileNotFoundException e)
			{
				System.err.println("Could not find the file: " + read_directory + read_name);
				e.printStackTrace();
			}
		formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
		formatter.setAddReturn(true);
		parser = new PDB_File_Parser();

	}


	/* ************************** METHODS ******************************************************************* */

	/**
	 * @return This method returns a vector of Atoms obtained from the PDB File Parser; Reads the whole PDB file.
	 */
	public Vector<Atom> Read_PDB()
		{
			try
				{
					atoms = parser.parse_PDB(pdb_reader, formatter);
					atoms_CA = parser.get_Alpha_Carbons();
					chain_ids_read = parser.get_chain_IDs_Read();
					residues_read = parser.get_Residues_Read();
					residue_id_pairs_read = parser.get_Residue_ID_Pairs();

				} catch (FileNotFoundException e)
				{
					System.err.println("'File Not Found' Exception Thrown: " + read_directory + read_name);
					e.printStackTrace();
					System.exit(0);
				} catch (NumberFormatException e)
				{
					System.err.println("'Number Format' Exception Thrown. The format of the PDB file is non-standard.");
					e.printStackTrace();
					System.exit(0);
				} catch (IOException e)
				{
					System.err.println("'IO' Exception Thrown. Could not read the file: " + read_directory + read_name);
					e.printStackTrace();
					System.exit(0);
				} catch (ArrayIndexOutOfBoundsException e)
				{
					System.err.println("'Array Index Out Of Bounds' Exception thrown for PDB file: " + read_directory + read_name);
					System.err.println("Check that all PDB files in the working directory have the EXACT same number of residues (alpha carbons).");
					e.printStackTrace();
					System.exit(0);
				} catch (Exception e)
				{
					System.err.println("Exception thrown for PDB file: " + read_directory + read_name);
					System.err.println("Carefully Check the FORMAT of the PDB file.");
					e.printStackTrace();
					System.exit(0);
				}
			return atoms;
		}
	
	/**
	 * @return A vector of Atoms obtained from the PDB File Parser; Reads a whole PDB file.
	 */
	public Vector<Atom> Read_PDB(String path_to_PDB)
		{
			path = path_to_PDB;
			pdb = new File(path);
			try
				{
					pdb_reader = new BufferedReader(new FileReader(pdb));
				} catch (FileNotFoundException e)
				{
					System.err.println("Could not find the file: " + path);
					e.printStackTrace();
				}
			formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
			formatter.setAddReturn(true);
			parser = new PDB_File_Parser();
			try
				{
					atoms = parser.parse_PDB(pdb_reader, formatter);

				} catch (FileNotFoundException e)
				{
					System.err.println("'File Not Found' Exception Thrown: " + path);
					e.printStackTrace();
					System.exit(0);
				} catch (NumberFormatException e)
				{
					System.err.println("'Number Format' Exception Thrown. The format of the PDB file is non-standard.");
					e.printStackTrace();
					System.exit(0);
				} catch (IOException e)
				{
					System.err.println("'IO' Exception Thrown. Could not read the file: " + path);
					e.printStackTrace();
					System.exit(0);
				} catch (ArrayIndexOutOfBoundsException e)
				{
					System.err.println("'Array Index Out Of Bounds' Exception thrown for PDB file: " + path);
					System.err.println("Check that all PDB files in the working directory have the EXACT same number of residues (alpha carbons).");
					e.printStackTrace();
					System.exit(0);
				} catch (Exception e)
				{
					System.err.println("Exception thrown for PDB file: " + path);
					System.err.println("Carefully Check the FORMAT of the PDB file.");
					e.printStackTrace();
					System.exit(0);
				}
			return atoms;
		}

	/**
	 * This method extracts a Subset of residues from Single Chain PDB files.
	 * 
	 * @param res
	 *            The list of residue numbers to read.
	 * @return A vector of Atoms obtained from the PDB File Parser; Reads the specified list of residues.
	 */
	public Vector<Atom> Read_PDB_Subset(List<Integer> res)
		{

			try
				{
					atoms = parser.parse_PDB(pdb_reader, formatter, res);

				} catch (FileNotFoundException e)
				{
					System.err.println("'File Not Found' Exception Thrown: " + read_directory + read_name);
					e.printStackTrace();
					System.exit(0);
				} catch (NumberFormatException e)
				{
					System.err.println("'Number Format' Exception Thrown. The format of the PDB file is non-standard.");
					e.printStackTrace();
					System.exit(0);
				} catch (IOException e)
				{
					System.err.println("'IO' Exception Thrown. Could not read the file: " + read_directory + read_name);
					e.printStackTrace();
					System.exit(0);
				} catch (ArrayIndexOutOfBoundsException e)
				{
					System.err.println("'Array Index Out Of Bounds' Exception thrown for PDB file: " + read_directory + read_name);
					System.err.println("Check that all PDB files in the working directory have the EXACT same number of residues (alpha carbons).");
					e.printStackTrace();
					System.exit(0);
				} catch (Exception e)
				{
					System.err.println("Exception thrown for PDB file: " + read_directory + read_name);
					System.err.println("Carefully Check the FORMAT of the PDB file.");
					e.printStackTrace();
					System.exit(0);
				}
			return atoms;
		}

	/**
	 * This method extracts a Subset of residues from Multi Chain PDB files.
	 * 
	 * @param chainIDs
	 *            The list of chain IDs to read
	 * @param res
	 *            The list of residue numbers to read
	 * @return A vector of Atoms obtained from the PDB File Parser; Reads the specified list of chain IDs and residues.
	 */
	public Vector<Atom> Read_PDB_Subset_Multi(List<String> chainIDs, List<Integer> res)
		{
			try
				{
					atoms = parser.parse_PDB(pdb_reader, formatter, chainIDs, res);
				} catch (FileNotFoundException e)
				{
					System.err.println("'File Not Found' Exception Thrown: " + read_directory + read_name);
					e.printStackTrace();
					System.exit(0);
				} catch (NumberFormatException e)
				{
					System.err.println("'Number Format' Exception Thrown. The format of the PDB file is non-standard.");
					e.printStackTrace();
					System.exit(0);
				} catch (IOException e)
				{
					System.err.println("'IO' Exception Thrown. Could not read the file: " + read_directory + read_name);
					e.printStackTrace();
					System.exit(0);
				} catch (ArrayIndexOutOfBoundsException e)
				{
					System.err.println("'Array Index Out Of Bounds' Exception thrown for PDB file: " + read_directory + read_name);
					System.err.println("Check that all PDB files in the working directory have the EXACT same number of residues (alpha carbons).");
					e.printStackTrace();
					System.exit(0);
				} catch (Exception e)
				{
					System.err.println("Exception thrown for PDB file: " + read_directory + read_name);
					System.err.println("Carefully Check the FORMAT of the PDB file.");
					e.printStackTrace();
					System.exit(0);
				}
			return atoms;
		}

	/**
	 * This method writes a set of Atoms to a PDB file
	 * 
	 * @param w_dir
	 *            The write directory
	 * @param w_name
	 *            The name of the PDB file
	 * @param e_atoms
	 *            The set of Atoms
	 */
	public void Write_PDB(String w_dir, String w_name, Vector<Atom> e_atoms)
		{

			write_directory = w_dir;
			write_name = w_name;
			edited_atoms = e_atoms;
			edited_pdb = new File(write_directory + write_name);
			try
				{
					pdb_writer = new BufferedWriter(new FileWriter(edited_pdb));
					parser.write_PDB(pdb_writer, edited_atoms, formatter);
					pdb_writer.close();
				} catch (IOException io)
				{
					System.err.println("IOException Thrown. Could not write the file: " + write_directory + write_name);
					io.printStackTrace();
				}
		}

	/* ************************** GETTERS ******************************************************************* */

	/**
	 * @return Returns the set of alpha carbons when reading the entire PDB file.
	 */
	public Vector<Atom> get_Atoms_CA()
		{
			return atoms_CA;
		}

	/**
	 * @return The chain IDs read
	 */
	public List<String> getChain_ids_read()
		{
			return chain_ids_read;
		}

	/**
	 * @return The residue numbers read
	 */
	public List<Integer> getResidues_read()
		{
			return residues_read;
		}

	/**
	 * @return The residue_id_pairs_read (Chain ID + Residue Number)
	 */
	public List<String> getResidue_id_pairs_read()
		{
			return residue_id_pairs_read;
		}
}
