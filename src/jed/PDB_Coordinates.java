package jed;

import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import Jama.Matrix;

/**
 * JED class PDB_Coordinates: Data structure for the Alpha Carbon Coordinates in a PDB file.
 * Note: PCA is only done for the Alpha Carbons in this version.
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

class PDB_Coordinates
{

	int number_of_alpha_carbons;
	String directory, pdb_file;
	double x, y, z;
	Matrix col_vector;
	List<Double> x_coords, y_coords, z_coords;
	Vector<Atom> atoms, atoms_CA;

	/* ******************************* CONSTRUCTOR **************************************************************************** */

	/**
	 * Constructor to create a PDB Coordinates data structure.
	 * 
	 * @param dir
	 *            The working directory
	 * @param pdb_file
	 *            The source PDB file
	 */
	PDB_Coordinates(String dir, String pdb_file)
		{

			this.directory = dir;
			this.pdb_file = pdb_file;

			x_coords = new ArrayList<Double>();
			y_coords = new ArrayList<Double>();
			z_coords = new ArrayList<Double>();

			atoms = new Vector<Atom>();
			atoms_CA = new Vector<Atom>();

			get_PDB_Coords(pdb_file);
		}

	/* ******************************** METHODS ************************************************************************** */

	/**
	 * Method to extract the Cartesian coordinates {x,y,z} from the specified PDB file.
	 * 
	 * @param pdb_file
	 *            The source PDB file
	 */
	private void get_PDB_Coords(String pdb_file)
		{
			PDB_IO pdb = new PDB_IO(directory, pdb_file);
			atoms = pdb.Read_PDB();
			atoms_CA = pdb.get_Atoms_CA();

			number_of_alpha_carbons = 0;

			for (Atom a : atoms_CA)
				{
					x = a.x;
					y = a.y;
					z = a.z;

					x_coords.add(x);
					y_coords.add(y);
					z_coords.add(z);

					number_of_alpha_carbons++;
				}

			col_vector = new Matrix(3 * number_of_alpha_carbons, 1);

			for (int i = 0; i < number_of_alpha_carbons; i++)
				{
					col_vector.set(i, 0, x_coords.get(i));
					col_vector.set((i + number_of_alpha_carbons), 0, y_coords.get(i));
					col_vector.set((i + (2 * number_of_alpha_carbons)), 0, z_coords.get(i));
				}
		}

	/* ********************************* GETTERS ********************************************************************************* */

	/**
	 * @return The column vector of the Cartesian coordinates {x,y,z} extracted from the PDB file.
	 */
	public Matrix get_column_vector()
		{

			return col_vector;

		}

	/**
	 * @return The number of alpha carbons in the PDB file (Equal to the Column Vector Row_dimension/3)
	 */
	public double number_of_alpha_carbons()
		{

			return number_of_alpha_carbons;

		}
}
