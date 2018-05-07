package jed;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import Jama.Matrix;

/**
 * JED class Remove_Outliers_by_Frame: Class for removing outliers based on conformation RMSD of frames prior to PCA. Note: The expected input is 2 matrices of eigenvectors
 * representing 2 equidimensional subspaces from a vector space. If you are not sure if you are working with orthonormal bases (subsets of an eigenspace), there is a
 * test_Orthogonal method. Orthonormal bases are expected. This means that the number of rows and columns must be the same. Copyright (C) 2012 Dr. Charles David
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Dr. Charles David
 * 
 */

public class Remove_Outliers_by_Frame
{

	int COLS, ROWS, trim, cols_trimmed, cols_removed, length;
	double percent;
	List<Double> conformation_rmsds;
	List<Integer> removed_conformations;
	Matrix coordinates, coordinates_trimmed_by_frame, means, sum_of_sq_dev, sigmas, z_scores;

	/* **************************************** CONSTRUCTOR *************************************************************************** */

	/**
	 * Constructor to handle outliers by removing frames based on conformation RMSD.
	 * 
	 * @param input_data
	 *            The coordinates matrix
	 * @param conf_rmsds
	 *            The list of conformation RMSDs for the trajectory
	 */
	Remove_Outliers_by_Frame(Matrix input_data, List<Double> conf_rmsds)
		{

			coordinates = input_data;
			conformation_rmsds = conf_rmsds;

			COLS = coordinates.getColumnDimension();
			ROWS = coordinates.getRowDimension();

			length = conformation_rmsds.size();
		}

	/* **************************************** METHODS *************************************************************************** */

	/**
	 * Method that determines which frames are outliers based on the conformation RMSDs and removes them.
	 */
	public void remove_Frames()
		{
			double[] crmsds = new double[length];
			int i = 0;
			for (double d : conformation_rmsds)
			{
				crmsds[i] = d;
				i++;
			}
			double mean = Descriptive_Stats.get_standard_deviation(crmsds);
			double ss = Descriptive_Stats.get_sum_of_squared_deviations(crmsds, mean);
			double std_dev = Descriptive_Stats.get_standard_deviation(crmsds, mean, ss);
			double[] z = Descriptive_Stats.get_Z_scores(crmsds, mean, std_dev);
			z_scores = new Matrix(z, length);

			trim = (int) (COLS * percent);
			cols_trimmed = COLS - (trim);
			coordinates_trimmed_by_frame = new Matrix(ROWS, cols_trimmed);

			removed_conformations = new ArrayList<>();

			List<Double> sorted_data = new ArrayList<>();
			sorted_data.addAll(conformation_rmsds);
			Collections.sort(sorted_data, Collections.reverseOrder());

			// remove the high rmsds
			for (int j = 0; j < trim; j++)
			{
				double rv = sorted_data.get(0);
				int conf_num = conformation_rmsds.indexOf(rv);
				removed_conformations.add(conf_num);
				sorted_data.remove(0);
			}
			// Sort the list of removed conformations
			Collections.sort(removed_conformations);

			// create the trimmed matrix
			for (int j = 0; j < cols_trimmed; j++)
			{
				double value = sorted_data.get(j);
				int index = conformation_rmsds.indexOf(value);
				Matrix col = coordinates.getMatrix(0, ROWS - 1, index, index);
				coordinates_trimmed_by_frame.setMatrix(0, ROWS - 1, j, j, col);
			}
		}

	/* **************************************** SETTERS *************************************************************************** */

	/**
	 * Method to set the percentage of frames to remove.
	 */
	public void set_percent(double cut_percent)
		{
			this.percent = cut_percent;
		}

	/* **************************************** GETTERS *************************************************************************** */

	/**
	 * @return The list of all frames removed
	 */
	public List<Integer> get_Removed_Conformation_List()
		{

			return removed_conformations;
		}

	/**
	 * @return The Z scores for each frame (conformation)
	 */
	public Matrix get_Conformation_Z_Scores()
		{

			return z_scores;
		}

	/**
	 * @return The number of frames (columns) removed
	 */
	public int get_cols_removed()
		{
			return cols_removed;
		}

	/**
	 * @return The coordinates_trimmed_by_frame
	 */
	public Matrix get_Coordinates_Trimmed()
		{
			return coordinates_trimmed_by_frame;
		}
}
