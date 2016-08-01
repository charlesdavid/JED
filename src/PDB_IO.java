package jed;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Vector;

public class PDB_IO
{
	public static String path;
	public static Vector<Atom> atoms;
	public static File pdb;
	public static BufferedReader pdb_reader;
	public static BufferedWriter pdb_writer;
	public static FortranFormat formatter;
	public static PDB_File_Parser parser;

	/**
	 * @return A vector of Atoms obtained from the PDB File Parser; Reads a whole PDB file.
	 */
	public static Vector<Atom> Read_PDB(String path_to_PDB)
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
	 * This method writes a set of Atoms to a PDB file
	 * 
	 * @param w_dir
	 *            The write directory
	 * @param w_name
	 *            The name of the PDB file
	 * @param e_atoms
	 *            The set of Atoms
	 */
	public static void Write_PDB(String write_path, Vector<Atom> set_of_atoms)
		{

			path = write_path;
			pdb = new File(path);
			atoms = set_of_atoms;
			try
				{
					pdb_writer = new BufferedWriter(new FileWriter(pdb));
					parser.write_PDB(pdb_writer, atoms, formatter);
					pdb_writer.close();
				} catch (IOException io)
				{
					System.err.println("IOException Thrown. Could not write the file: " + path);
					io.printStackTrace();
				}
		}
}
