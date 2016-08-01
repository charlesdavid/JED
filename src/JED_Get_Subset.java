package jed;

import java.util.ArrayList;

import Jama.Matrix;

/**
 * JED class JED_Get_Subset: Constructs the subset of coordinates as specified in a residue list.
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
 * @author Dr. Charles David
 * 
 */

public class JED_Get_Subset
{

	static String directory, out_dir, description;
	static int number_of_residues, number_of_residues_SS, number_of_conformations;
	static int ROWS, ROWS_SS, COLS, COLS_SS;
	static double trace;
	static int[] residue_list;
	static ArrayList<Integer> res_list;
	static Matrix original_PDB_coordinates, subset_PDB_coordinates;

	// ***************************************** CONSTRUCTORS ***************************************************//

	/**
	 * Constructor for extracting a subset of residue coordinates from the original PDB coordinates matrix
	 * 
	 * @param dir
	 *            The working directory
	 * @param des
	 *            The job description
	 * @param data
	 *            The original PDB coordinates matrix
	 * @param res_list
	 *            The list of residues that defines the subset
	 */
	JED_Get_Subset(String dir, String des, Matrix data, int[] res_list)
		{

			directory = dir;
			description = des;
			out_dir = directory + "JED_RESULTS_" + description + "/";
			original_PDB_coordinates = data;
			residue_list = res_list;
			ROWS = original_PDB_coordinates.getRowDimension();
			COLS = original_PDB_coordinates.getColumnDimension();
			number_of_residues = ROWS / 3;
			number_of_residues_SS = residue_list.length;
			ROWS_SS = number_of_residues_SS * 3;
			COLS_SS = COLS;
			number_of_conformations = COLS;
			subset_PDB_coordinates = new Matrix(ROWS_SS, COLS_SS);
		}

	// ***************************************** METHODS ***************************************************//

	/**
	 * Calculates the subset matrix of PDB coordinates
	 * 
	 * @return Subset Coordinates Matrix
	 */
	Matrix get_Subset_Coords()
		{

			for (int i = 0; i < number_of_residues_SS; i++) // array indices must reference ZERO as first element (Java Array Numbering)
				{
					for (int j = 0; j < number_of_conformations; j++)
						{
							double element_X = original_PDB_coordinates.get((residue_list[i]), j);
							double element_Y = original_PDB_coordinates.get((residue_list[i] + (number_of_residues)), j);
							double element_Z = original_PDB_coordinates.get((residue_list[i] + (2 * number_of_residues)), j);

							subset_PDB_coordinates.set(i, j, element_X);
							subset_PDB_coordinates.set((i + number_of_residues_SS), j, element_Y);
							subset_PDB_coordinates.set((i + (2 * number_of_residues_SS)), j, element_Z);
						}
				}
			String name = "ss_" + number_of_residues_SS + "_PDB_coordinates.txt";
			Matrix_IO.write_Matrix(subset_PDB_coordinates, out_dir, name, 9, 3);
			return subset_PDB_coordinates;
		}

	/**
	 * Sets the output directory
	 * 
	 * @param dir
	 *            The directory for output
	 */
	public void set_Output_Directory(String dir)
		{
			out_dir = dir;
		}
}
