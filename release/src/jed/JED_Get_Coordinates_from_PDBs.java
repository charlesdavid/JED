package jed;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import Jama.Matrix;

/**
 * JED class JED_Get_Coordinate_from_PDBs: Reads all PDB files in a specified directory to create the Coordinates Matrix.
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
 * along with this program. If not, see <http://www.gnu.org/license>
 * 
 * @ author Dr. Charles David
 */

public class JED_Get_Coordinates_from_PDBs
{

	static int number_of_residues, number_of_conformations, ROWS, COLS;
	static String directory, out_dir, description, pdb_data_file;
	static Matrix PDB_coordinates, col_vector;
	static File pdb_file, orig_PDB_coords;
	static FilenameFilter filter;
	static String[] ls;
	static PDB_Coordinates data;
	static List<String> file_names;

	/**
	 * Constructor to get the matrix of PDB coordinates by reading all the PDB files in the specified directory.
	 * 
	 * @param dir
	 *            The directory from which to read the PDB files
	 * @param des
	 *            The job description
	 */
	@SuppressWarnings("unused")
	public JED_Get_Coordinates_from_PDBs(String dir, String des)
		{

			directory = dir;
			description = des;
			out_dir = directory + "JED_RESULTS_" + description + "/";
			boolean exist = new File(out_dir).exists();
			if (!exist)
				{
					boolean success = (new File(out_dir)).mkdirs();
				}
			pdb_file = new File(directory);
			filter = new PDB_Filter();
			ls = pdb_file.list(filter);
			Arrays.sort(ls);
			data = new PDB_Coordinates(directory, ls[0]);
			col_vector = data.get_column_vector();
			ROWS = col_vector.getRowDimension();
			number_of_residues = ROWS / 3;
			COLS = ls.length;
			number_of_conformations = ls.length;
			file_names = new ArrayList<String>();
			PDB_coordinates = new Matrix(ROWS, COLS);
		}

	/**
	 * Returns the matrix of PDB coordinates
	 * 
	 * @return The original PDB Coordinates Matrix
	 */
	public Matrix get_Orig_PDB_Coords()
		{

			for (int i = 0; ls != null && i < ls.length;)
				for (i = 0; i < ls.length; i++)
					{
						data = new PDB_Coordinates(directory, ls[i]);
						file_names.add(ls[i]);
						Matrix cv = data.get_column_vector();
						PDB_coordinates.setMatrix(0, ROWS - 1, i, i, cv);
						data = null;
						if ((i % 10 == 0)) System.gc();
					}
			String path = out_dir + "original_PDB_Coordinates.txt";
			Matrix_IO.write_Matrix(PDB_coordinates, path, 9, 3);
			return PDB_coordinates;
		}

	/**
	 * Returns the list of all the PDB files read, in the order in which they were read
	 * 
	 * @return The list of PDB file names
	 */
	public List<String> get_PDB_names()
		{

			return file_names;
		}
}
