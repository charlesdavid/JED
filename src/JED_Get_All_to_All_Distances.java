package jed;

import java.io.File;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.List;

import Jama.Matrix;

/**
 * JED class Get_All_to_All_Distances: Gets the n(n-1)/2 matrix of distances for n residues.
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

public class JED_Get_All_to_All_Distances
{

	static String directory, description, out_dir;
	static int ROWS, COLS, number_of_residues_dist;
	static Matrix subset_coords, distances, distances_T;
	static NumberFormat nf;
	static RoundingMode rm;
	static List<String> chain_IDs;
	static List<Integer> chain_offsets;

	/**
	 * Constructor for getting the All-to-All Distances
	 * 
	 * @param data
	 *            The Subset Coordinates Matrix
	 * @param dir
	 *            The working directory
	 * @param desc
	 *            The job description
	 */
	@SuppressWarnings("unused")
	public JED_Get_All_to_All_Distances(Matrix data, String dir, String desc)
		{

			subset_coords = data;
			ROWS = subset_coords.getRowDimension();
			number_of_residues_dist = (ROWS / 3);
			directory = dir;
			description = desc;
			out_dir = directory + "JED_RESULTS_" + description + "/dPCA/";
			boolean exist = new File(out_dir).exists();
			if (!exist)
				{
					boolean success = (new File(out_dir)).mkdirs();
				}

			nf = NumberFormat.getInstance();
			rm = RoundingMode.HALF_UP;
			nf.setRoundingMode(rm);
			nf.setMaximumFractionDigits(3);
			nf.setMinimumFractionDigits(3);
		}

	/* *************************************************************************************** */

	/**
	 * Gets the Matrix of the All to All Distances, and writes the transpose of the matrix to file.
	 * 
	 * @return
	 */
	public Matrix get_All_To_All_Distances()
		{

			distances = do_distances();
			distances_T = distances.transpose();
			String path = out_dir + "ss_" + number_of_residues_dist + "_all_to_all_distances.txt";
			Matrix_IO.write_Matrix(distances_T, path, 9, 3);
			return distances;
		}

	private Matrix do_distances()
		{

			int COLS = subset_coords.getColumnDimension();
			int ROWS_D = ((number_of_residues_dist * (number_of_residues_dist - 1)) / 2);
			Matrix dists = new Matrix(ROWS_D, COLS);
			for (int observation = 0; observation < COLS; observation++)
				{
					int count = 0;
					for (int i = 0; i < number_of_residues_dist; i++)
						{
							double x = subset_coords.get(i, observation);
							double y = subset_coords.get(i + number_of_residues_dist, observation);
							double z = subset_coords.get(i + 2 * number_of_residues_dist, observation);

							for (int j = (i + 1); j < number_of_residues_dist; j++, count++)
								{

									double xr = subset_coords.get((j), observation);
									double yr = subset_coords.get((j) + number_of_residues_dist, observation);
									double zr = subset_coords.get((j) + 2 * number_of_residues_dist, observation);

									double xxr = Math.pow((x - xr), 2);
									double yyr = Math.pow((y - yr), 2);
									double zzr = Math.pow((z - zr), 2);

									double sum_xyz = xxr + yyr + zzr;
									double d = Math.sqrt(sum_xyz);
									dists.set(count, observation, d);
								}
						}
				}
			return dists;
		}
}
