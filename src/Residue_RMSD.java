package jed;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;

/**
 * JED class Residue_RMSD: This class computes the Residue RMSDs (RMSFs) from an entire trajectory.
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

public class Residue_RMSD
{

	int COLS, ROWS, number_of_alpha_carbons;
	Matrix coordinates, means, sum_of_sq_dev, sigmas, z_scores;

	/* ************************************** CONSTRUCTORS ******************************************************************************** */

	Residue_RMSD(Matrix input)
		{

			coordinates = input;
			COLS = coordinates.getColumnDimension();
			ROWS = coordinates.getRowDimension();
			number_of_alpha_carbons = (ROWS / 3);
			means = new Matrix(ROWS, 1);
			sum_of_sq_dev = new Matrix(ROWS, 1);
			sigmas = new Matrix(ROWS, 1);
			z_scores = new Matrix(ROWS, COLS);
		}

	/* ************************************** METHODS ******************************************************************************** */

	private void get_residue_stats()
		{

			for (int i = 0; i < ROWS; i++)
				{
					double[] row = coordinates.getMatrix(i, i, 0, COLS - 1).getRowPackedCopy();
					double mean = Descriptive_Stats.get_mean(row);
					double ssdevs = Descriptive_Stats.get_sum_of_squared_deviations(row, mean);
					double sigma = Descriptive_Stats.get_standard_deviation(row, mean, ssdevs);
					means.set(i, 0, mean);
					sum_of_sq_dev.set(i, 0, ssdevs);
					sigmas.set(i, 0, sigma);
					double[] z_s = Descriptive_Stats.get_Z_scores(row, mean, sigma);
					Matrix z = new Matrix(z_s, 1);
					z_scores.setMatrix(i, i, 0, COLS - 1, z);
				}
		}

	private Matrix get_residue_rmsds_matrix()
		{

			get_residue_stats();

			Matrix res_rmsds = new Matrix(number_of_alpha_carbons, 1);
			for (int i = 0; i < number_of_alpha_carbons; i++)
				{
					double x = sum_of_sq_dev.get(i, 0);
					double y = sum_of_sq_dev.get(i + number_of_alpha_carbons, 0);
					double z = sum_of_sq_dev.get(i + 2 * number_of_alpha_carbons, 0);
					double sum = ((x + y + z) / COLS);
					double rmsd = Math.sqrt(sum);
					res_rmsds.set(i, 0, rmsd);
				}
			return res_rmsds;
		}

	/**
	 * @return The list of residue rmsds (RMSF) based on the entire trajectory.
	 */
	public List<Double> get_residue_rmsd()
		{

			Matrix rmsds = get_residue_rmsds_matrix();

			ArrayList<Double> residue_rmsds = new ArrayList<Double>();
			double[] rrmsds = rmsds.getColumnPackedCopy();
			for (double d : rrmsds)
				{
					residue_rmsds.add(d);
				}
			return residue_rmsds;
		}

	/**
	 * @return The variable Z scores.
	 */
	public Matrix get_z_scores()
		{
			return z_scores;
		}
}
