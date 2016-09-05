package jed;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

/**
 * JED class PDB_File_Parser: Parser class to read and write PDb files using Fortran Format.
 * The exceptions thrown by this class are caught by the JED_Read_PDBs class.
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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Dr. Charles David
 */

class PDB_File_Parser
{

	String line, CA = "CA";
	StringBuilder sb1, sb2;
	List<Integer> residue_list, residues_read;
	List<String> lines, chainID_list, chain_IDs_read, residue_ID_pairs_read;
	Vector<Atom> atoms, atoms_CA;

	/* **************************** CONSTRUCTORS **************************************** */

	/**
	 * Constructs a PDB File Parser object, initialized all array lists.
	 */
	public PDB_File_Parser()
		{
			chain_IDs_read = new ArrayList<String>();
			chainID_list = new ArrayList<String>();
			residue_list = new ArrayList<Integer>();
			residues_read = new ArrayList<Integer>();
			residue_ID_pairs_read = new ArrayList<String>();
			lines = new ArrayList<String>();
		}

	/* ******************************* METHODS ******************************************* */

	/**
	 * Reads a PDB file using the specified reader and formatter.
	 * 
	 * @param br
	 *            The buffered reader
	 * @param formatter
	 *            The Fortran Format formatter
	 * @return The Vector of Atoms
	 * @throws IOException
	 *             If the PDB file can not be read
	 */
	public Vector<Atom> parse_PDB(BufferedReader br, FortranFormat formatter) throws IOException
		{
			atoms = new Vector<Atom>();
			atoms_CA = new Vector<Atom>();
			sb1 = new StringBuilder();
			sb2 = new StringBuilder();

			while ((line = br.readLine()) != null)
				lines.add(line);
			br.close();

			for (String file_line : lines)
				{
					if (file_line.startsWith("SEQRES")) sb1.append(file_line + "\n");

					if (file_line.startsWith("ATOM"))
						{
							Vector<Object> objects = formatter.parse(file_line);
							Atom a = new Atom();
							a.header = (String) objects.get(0);
							a.atom_number = (Integer) objects.get(1);
							a.symbol = (String) objects.get(2);
							a.res_type = (String) objects.get(4);
							a.chainID = (String) objects.get(5);
							a.res_number = (Integer) objects.get(6);
							a.x = (Double) objects.get(8);
							a.y = (Double) objects.get(9);
							a.z = (Double) objects.get(10);
							if (objects.get(11) != null) a.occupancy = (Double) objects.get(11);
							a.b_factor = (Double) objects.get(12);
							a.code = (String) objects.get(13);

							if (a.symbol.equals(CA))
								{
									chain_IDs_read.add(a.chainID);
									residues_read.add(a.res_number);
									atoms_CA.add(a);
								}
							atoms.add(a);
						}

					if (file_line.startsWith("TER")) sb2.append(file_line + "\n");

					if (file_line.startsWith("CONECT")) sb2.append(file_line + "\n");
				}
			return atoms;
		}

	/**
	 * Reads a PDB file using the specified reader, formatter, and list of residue numbers
	 * 
	 * @param br
	 *            The buffered reader
	 * @param formatter
	 *            The Fortran Format formatter
	 * @param residues
	 *            The list of residue numbers to read
	 * @return The Vector of Atoms
	 * @throws IOException
	 *             If the PDB file can not be read
	 */
	public Vector<Atom> parse_PDB(BufferedReader br, FortranFormat formatter, List<Integer> residues) throws IOException
		{
			atoms = new Vector<Atom>();
			sb1 = new StringBuilder();
			sb2 = new StringBuilder();

			residue_list = residues;
			int num_of_residues = residue_list.size();

			while ((line = br.readLine()) != null)
				lines.add(line);
			br.close();

			for (String file_line : lines)
				{
					if (file_line.startsWith("SEQRES")) sb1.append(file_line + "\n");

					if (file_line.startsWith("ATOM"))
						{
							Vector<Object> objects = formatter.parse(file_line);
							Atom a = new Atom();
							a.header = (String) objects.get(0);
							a.atom_number = (Integer) objects.get(1);
							a.symbol = (String) objects.get(2);
							a.res_type = (String) objects.get(4);
							a.chainID = (String) objects.get(5);
							a.res_number = (Integer) objects.get(6);
							a.x = (Double) objects.get(8);
							a.y = (Double) objects.get(9);
							a.z = (Double) objects.get(10);
							if (objects.get(11) != null) a.occupancy = (Double) objects.get(11);
							a.b_factor = (Double) objects.get(12);
							a.code = (String) objects.get(13);
							int key = a.res_number;
							for (int res_indx = 0; res_indx < num_of_residues; res_indx++)
								{
									if (residue_list.get(res_indx).equals(key)) atoms.add(a);
								}
						}
					if (file_line.startsWith("TER")) sb2.append(file_line + "\n");

					if (file_line.startsWith("CONECT")) sb2.append(file_line + "\n");
				}
			return atoms;
		}

	/**
	 * Reads a PDB file using the specified reader, formatter, list of chain IDs, and list of residue numbers
	 * 
	 * @param br
	 *            The buffered reader
	 * @param formatter
	 *            The Fortran Format formatter
	 * @param chain_ids
	 *            The list of residue chain IDs to read
	 * @param res
	 *            The list of residue numbers to read
	 * @return The Vector of Atoms
	 * @throws IOException
	 *             If the PDB file can not be read
	 */
	public Vector<Atom> parse_PDB(BufferedReader br, FortranFormat formatter, List<String> chain_ids, List<Integer> res) throws IOException
		{
			atoms = new Vector<Atom>();
			sb1 = new StringBuilder();
			sb2 = new StringBuilder();

			chainID_list = chain_ids;
			residue_list = res;
			int num_of_residues = residue_list.size();

			while ((line = br.readLine()) != null)
				lines.add(line);
			br.close();

			for (String file_line : lines)
				{
					if (file_line.startsWith("SEQRES")) sb1.append(file_line + "\n");

					if (file_line.startsWith("ATOM"))
						{
							Vector<Object> objects = formatter.parse(file_line);
							Atom a = new Atom();
							a.header = (String) objects.get(0);
							a.atom_number = (Integer) objects.get(1);
							a.symbol = (String) objects.get(2);
							a.res_type = (String) objects.get(4);
							a.chainID = (String) objects.get(5);
							a.res_number = (Integer) objects.get(6);
							a.x = (Double) objects.get(8);
							a.y = (Double) objects.get(9);
							a.z = (Double) objects.get(10);
							if (objects.get(11) != null) a.occupancy = (Double) objects.get(11);
							a.b_factor = (Double) objects.get(12);
							a.code = (String) objects.get(13);

							String ID = a.chainID;
							int key = a.res_number;

							for (int res_indx = 0; res_indx < num_of_residues; res_indx++)
								{
									if (chainID_list.get(res_indx).equals(ID) && residue_list.get(res_indx).equals(key)) atoms.add(a);
								}
						}

					if (file_line.startsWith("TER")) sb2.append(file_line + "\n");

					if (file_line.startsWith("CONECT")) sb2.append(file_line + "\n");
				}
			return atoms;
		}

	/* ------------------------------------------------------------------------------------ */

	/**
	 * Writes a PDB file using the specified writer, set of atoms, and formatter.
	 * 
	 * @param writer
	 *            The Buffered FileWriter
	 * @param atoms
	 *            The Vector of Atoms to write
	 * @param formatter
	 *            The Fortran Formatter
	 * @throws IOException
	 *             If the PDB file can not be written
	 */
	public void write_PDB(BufferedWriter writer, Vector<Atom> atoms, FortranFormat formatter) throws IOException
		{
			Vector<Object> objects = new Vector<Object>(15);
			if (sb1 != null) writer.write(sb1.toString());
			for (Atom a : atoms)
				{
					objects.clear();
					objects.add(a.header + (a.header.length() == 1 ? "   " : "  "));
					objects.add(a.atom_number);
					objects.add(a.symbol + (a.symbol.length() == 1 ? "   " : "  "));
					objects.add("");
					objects.add(a.res_type);
					objects.add(a.chainID);
					objects.add(a.res_number);
					objects.add(null);
					objects.add(a.x);
					objects.add(a.y);
					objects.add(a.z);
					objects.add(a.occupancy);
					objects.add(a.b_factor);
					objects.add(a.code);
					writer.append(formatter.format(objects));
				}
			if (sb2 != null) writer.write(sb2.toString());
			writer.write("END");
			writer.close();
		}

	/* ******************************* GETTERS ******************************************* */

	/**
	 * @return The list of residue ID pairs: (Chain_ID + Residue_Number)
	 */
	public List<String> get_Residue_ID_Pairs()
		{
			for (int i = 0; i < residues_read.size(); i++)
				{
					String ID_Pair = chain_IDs_read.get(i) + residues_read.get(i);
					residue_ID_pairs_read.add(ID_Pair);
				}
			return residue_ID_pairs_read;
		}

	/**
	 * @return The list of chain_IDs that the parser read.
	 */
	public List<String> get_chain_IDs_Read()
		{

			return chain_IDs_read;
		}

	/**
	 * @return The list of residue numbers that the parser read
	 */
	public List<Integer> get_Residues_Read()
		{

			return residues_read;
		}

	/**
	 * @return The vector of alpha carbons that were read by the parser
	 */
	public Vector<Atom> get_Alpha_Carbons()
		{
			return atoms_CA;
		}
}
