package jed;

import java.io.File;

import Jama.Matrix;

/**
 * JED class JED_Get_Distances_for_Residue_Pairs: Constructs the matrix of distances for the Residue Pairs subset.
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
 */

public class JED_Get_Distances_for_Residue_Pairs
	{

		String directory, description, out_dir;
		int ROWS, ROWS_dp, COLS, number_of_residues, number_of_residues_pairs;
		Matrix X_vectors, ref_distances, distances;
		int[] residue_list1, residue_list2;
		boolean exist;

		/* ***************************************** CONSTRUCTOR ********************************************************************* */

		/**
		 * Constructor to get the distances for the specified residue pairs
		 * 
		 * @param data
		 *            The original PDB coordinates matrix
		 * @param residues1
		 *            The residue list specifying all the first elements of the pairs
		 * @param residues2
		 *            The residue list specifying all the second elements of the pairs
		 * @param dir
		 *            The working directory
		 * @param desc
		 *            The job description
		 */
		public JED_Get_Distances_for_Residue_Pairs(Matrix data, int[] residues1, int[] residues2, String dir, String desc)
			{

				X_vectors = data;
				ROWS = X_vectors.getRowDimension();
				COLS = X_vectors.getColumnDimension();
				number_of_residues = (ROWS / 3);
				directory = dir;
				description = desc;
				out_dir = directory + "JED_RESULTS_" + description + "/dpPCA/";
				exist = new File(out_dir).exists();
				if (!exist)
					{
						(new File(out_dir)).mkdirs();
					}

				residue_list1 = residues1;
				residue_list2 = residues2;
				number_of_residues_pairs = residue_list1.length;
				ROWS_dp = residue_list1.length;
			}

		/* ******************************************* METHODS ********************************************************************* */

		/**
		 * Calculates the distances matrix for the specified residue pairs in the reference structure and writes it to file
		 * 
		 * @return The reference residue pair distances matrix
		 */
		public Matrix Get_Ref_Distances()
			{

				ref_distances = new Matrix(ROWS_dp, 1);

				for (int i = 0; i < number_of_residues_pairs; i++)
					{
						int residue1 = residue_list1[i]; // Because residue lists start with 1 while Java arrays start with 0
						int residue2 = residue_list2[i];

						double x = X_vectors.get(residue1, 0);
						double y = X_vectors.get(residue1 + number_of_residues, 0);
						double z = X_vectors.get(residue1 + 2 * number_of_residues, 0);

						double xr = X_vectors.get(residue2, 0);
						double yr = X_vectors.get(residue2 + number_of_residues, 0);
						double zr = X_vectors.get(residue2 + 2 * number_of_residues, 0);

						double xxr = Math.pow((x - xr), 2);
						double yyr = Math.pow((y - yr), 2);
						double zzr = Math.pow((z - zr), 2);

						double sum_xyz = xxr + yyr + zzr;
						double d = Math.sqrt(sum_xyz);
						ref_distances.set(i, 0, d);
					}

				String path = out_dir + "ss_" + number_of_residues_pairs + "_Reference_Residue_Pairs_Distances.txt";
				Matrix_IO.write_Matrix(ref_distances.transpose(), path, 9, 3);
				return ref_distances;
			}

		/**
		 * Calculates the distances matrix for all the specified residue pairs and writes it to file
		 * 
		 * @return The residue pair distances matrix
		 */
		public Matrix get_Distances()
			{

				distances = new Matrix(ROWS_dp, COLS);
				for (int observation = 0; observation < COLS; observation++)
					{
						for (int i = 0; i < number_of_residues_pairs; i++)
							{
								int residue1 = residue_list1[i]; // Because residue lists start with 1 while Java arrays start with 0
								int residue2 = residue_list2[i];

								double x = X_vectors.get(residue1, observation);
								double y = X_vectors.get(residue1 + number_of_residues, observation);
								double z = X_vectors.get(residue1 + 2 * number_of_residues, observation);

								double xr = X_vectors.get(residue2, observation);
								double yr = X_vectors.get(residue2 + number_of_residues, observation);
								double zr = X_vectors.get(residue2 + 2 * number_of_residues, observation);

								double xxr = Math.pow((x - xr), 2);
								double yyr = Math.pow((y - yr), 2);
								double zzr = Math.pow((z - zr), 2);

								double sum_xyz = xxr + yyr + zzr;
								double d = Math.sqrt(sum_xyz);
								distances.set(i, observation, d);
							}
					}
				String path = out_dir + "ss_" + number_of_residues_pairs + "_Residue_Pairs_Distances.txt";
				Matrix_IO.write_Matrix(distances.transpose(), path, 9, 3);
				return distances;
			}
	}
