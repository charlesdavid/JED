package jed;

import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;

/**
 * JED class JED_Get_Transformed_Coordinates: Optimally aligns all frames of a trajectory to a reference frame using Quaternion operations. Copyright (C) 2012 Dr. Charles David
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/license>
 * 
 * @author Dr. Charles David
 */

public class JED_Get_Transformed_Coordinates
{

	static String directory, out_dir, description;
	static double z_cut_off, percent_cut;
	static int ROWS, COLS, number_of_residues;
	static List<Double> original_conformation_rmsds, transformed_conformation_rmsds, transformed_residue_rmsd_list, Z_Scores;
	static List<Integer> conformation_outliers;
	// static double[] var_means;
	static Matrix subset_REF_PDB_coordinates, transformed_subset_REF_PDB_coordinates, subset_PDB_coordinates, transformed_subset_PDB_coordinates, trimmed_PDB_coordinates_cols,
			adjusted_PDB_coordinates_rows, var_Z_scores, conf_Z_scores;
	static NumberFormat nf;
	static RoundingMode rm;
	static JED_Transform_Coords tf_coords;

	/* ************************************** CONSTRUCTORS ******************************************************************************** */

	/**
	 * Constructor for transforming the original coordinates using quaternion operations and calculating the conformation and residue RMSDs. If specified, the coordinates outliers
	 * will be handled This class is used for Cartesian analysis, but not the distance analyses which use internal coordinates (distances)
	 * 
	 * @param data
	 *            The original PDB coordinates
	 * @param ref_pdb
	 *            The reference PDB structure to use for the transformation
	 * @param dir
	 *            The working directory
	 * @param des
	 *            The job description
	 */
	JED_Get_Transformed_Coordinates(Matrix data, Matrix ref_coords, String dir, String des)
		{

			nf = NumberFormat.getInstance();
			rm = RoundingMode.HALF_UP;
			nf.setRoundingMode(rm);
			nf.setMaximumFractionDigits(3);
			nf.setMinimumFractionDigits(3);

			subset_PDB_coordinates = data;
			subset_REF_PDB_coordinates = ref_coords;
			directory = dir;
			description = des;

			out_dir = directory + "JED_RESULTS_" + description + "/";
			transformed_residue_rmsd_list = new ArrayList<>();
			original_conformation_rmsds = new ArrayList<>();
			transformed_conformation_rmsds = new ArrayList<>();

			ROWS = subset_PDB_coordinates.getRowDimension();
			COLS = subset_PDB_coordinates.getColumnDimension();
			number_of_residues = (ROWS / 3);
			transformed_subset_PDB_coordinates = new Matrix(ROWS, COLS);

			tf_coords = new JED_Transform_Coords(subset_REF_PDB_coordinates, subset_REF_PDB_coordinates);
			transformed_subset_REF_PDB_coordinates = tf_coords.get_transformed_coords();
		}

	/* ************************************** SETTERS ******************************************************************************** */

	/**
	 * Sets the Z cutoff for identifying outliers
	 * 
	 * @param z
	 *            The Z-cutoff
	 */
	public void set_z_cutoff(double z)
		{

			z_cut_off = z;
		}

	/**
	 * Sets the percent of frames to treat as outliers
	 * 
	 * @param percent
	 *            The percent of frames to remove
	 */
	public void set_percent_cutoff(double percent)
		{

			percent_cut = percent;
		}

	/**
	 * Sets the output directory
	 * 
	 * @param dir
	 *            The directory to use for output
	 */
	public void set_Output_Directory(String dir)
		{

			out_dir = dir;
		}

	/* ************************************** METHODS ******************************************************************************** */

	/**
	 * @return The Transformed coordinates for the specified subset
	 */
	public Matrix get_SS_Transformed_coords()
		{
			for (int z = 0; z < COLS; z++)
			{
				Matrix fc = subset_PDB_coordinates.getMatrix(0, ROWS - 1, z, z);
				tf_coords = new JED_Transform_Coords(transformed_subset_REF_PDB_coordinates, fc);
				Matrix transformed_col_vector = tf_coords.get_transformed_coords();
				transformed_subset_PDB_coordinates.setMatrix(0, ROWS - 1, z, z, transformed_col_vector);
				tf_coords = null;
				if (z % 10 == 0) System.gc();
			}
			String path = out_dir + "ss_" + number_of_residues + "_transformed_PDB_coordinates.txt";
			Matrix_IO.write_Matrix(transformed_subset_PDB_coordinates, path, 9, 3);
			path = out_dir + "ss_" + number_of_residues + "_original_PDB_Reference_coordinates.txt";
			Matrix_IO.write_Matrix(subset_REF_PDB_coordinates, path, 9, 3);
			path = out_dir + "ss_" + number_of_residues + "_transformed_PDB_Reference_coordinates.txt";
			Matrix_IO.write_Matrix(transformed_subset_REF_PDB_coordinates, path, 9, 3);

			return transformed_subset_PDB_coordinates;
		}

	/**
	 * @return The transformed conformation RMSDs
	 */
	public List<Double> get_SS_Conformation_RMSDs()
		{
			for (int z = 0; z < COLS; z++)
			{
				Matrix fc_O = subset_PDB_coordinates.getMatrix(0, ROWS - 1, z, z);
				Matrix fc_T = transformed_subset_PDB_coordinates.getMatrix(0, ROWS - 1, z, z);

				JED_Get_RMSD gRMSD_O = new JED_Get_RMSD(subset_REF_PDB_coordinates, fc_O);
				JED_Get_RMSD gRMSD_T = new JED_Get_RMSD(transformed_subset_REF_PDB_coordinates, fc_T);

				double rmsd_O = gRMSD_O.get_RMSD();
				original_conformation_rmsds.add(rmsd_O);
				double rmsd_T = gRMSD_T.get_RMSD();
				transformed_conformation_rmsds.add(rmsd_T);
			}
			String path = out_dir + "ss_" + number_of_residues + "_original_conformation_rmsds.txt";
			List_IO.write_Double_List(original_conformation_rmsds, path, 3);
			path = out_dir + "ss_" + number_of_residues + "_transformed_conformation_rmsds.txt";
			List_IO.write_Double_List(transformed_conformation_rmsds, path, 3);

			return transformed_conformation_rmsds;
		}

	/**
	 * @return The PDB coordinates matrix with the outlier frames removed
	 */
	public Matrix get_SS_transformed_coordinates_trimmed_COLS()
		{
			if (percent_cut > 0)
			{
				Remove_Outliers_by_Frame tdc = new Remove_Outliers_by_Frame(transformed_subset_PDB_coordinates, transformed_conformation_rmsds);
				tdc.set_percent(percent_cut);
				tdc.remove_Frames();
				trimmed_PDB_coordinates_cols = tdc.get_Coordinates_Trimmed();
				double percentage = (int) (percent_cut * 100);
				String path = out_dir + "ss_" + number_of_residues + "_trimmed_" + percentage + "_percent_PDB_coordinates_COLS.txt";
				Matrix_IO.write_Matrix(trimmed_PDB_coordinates_cols, path, 9, 3);
				conformation_outliers = tdc.get_Removed_Conformation_List();
				path = out_dir + "ss_" + number_of_residues + "_Removed_Conformation_Outliers.txt";
				List_IO.write_Integer_List(conformation_outliers, path);
				conf_Z_scores = tdc.get_Conformation_Z_Scores();
				path = out_dir + "ss_" + number_of_residues + "_Conformation_Z_Scores.txt";
				Matrix_IO.write_Matrix(conf_Z_scores, path, 6, 1);
			} else
			{
				trimmed_PDB_coordinates_cols = transformed_subset_PDB_coordinates;
			}
			return trimmed_PDB_coordinates_cols;
		}

	/**
	 * @return The PDB coordinates matrix with the outliers adjusted to their mean values.
	 */
	public Matrix get_SS_transformed_coordinates_adjusted_ROWS()
		{

			if (z_cut_off > 0)
			{
				Adjust_Outliers_by_Z_Score adr = new Adjust_Outliers_by_Z_Score(trimmed_PDB_coordinates_cols);
				adr.set_Z_threshold(z_cut_off);
				adr.adjust_row_data();
				adjusted_PDB_coordinates_rows = adr.get_coorinates_adjusted();
				String path = out_dir + "ss_" + number_of_residues + "_Z_threshold_" + z_cut_off + "_adjusted_PDB_coordinates_ROWS.txt";
				Matrix_IO.write_Matrix(adjusted_PDB_coordinates_rows, path, 9, 3);
				Matrix var_counts = adr.get_counts();
				path = out_dir + "ss_" + number_of_residues + "_adjustments_per_variable.txt";
				Matrix_IO.write_Matrix(var_counts, path, 6, 0);
			} else
			{
				adjusted_PDB_coordinates_rows = trimmed_PDB_coordinates_cols;
			}
			return adjusted_PDB_coordinates_rows;
		}

	/**
	 * @return The Residue RMSDs (RMSFs)
	 */
	public List<Double> get_SS_Residue_RMSDs()
		{

			Residue_RMSD rrmsd = new Residue_RMSD(transformed_subset_PDB_coordinates);
			transformed_residue_rmsd_list = rrmsd.get_residue_rmsd();
			var_Z_scores = rrmsd.get_z_scores();

			String path = out_dir + "ss_" + number_of_residues + "_residue_rmsd.txt";
			List_IO.write_Double_List(transformed_residue_rmsd_list, path, 3);

			path = out_dir + "ss_" + number_of_residues + "_Variable_Z_Scores.txt";
			Matrix_IO.write_Matrix(var_Z_scores, path, 6, 1);

			return transformed_residue_rmsd_list;
		}

	/* ************************************** GETTERS ******************************************************************************** */

	/**
	 * @return the original_reference_coordinates
	 */
	public Matrix get_Original_reference_coordinates()
		{
			return subset_REF_PDB_coordinates;
		}

	/**
	 * @return the transformed_reference_coordinates
	 */
	public Matrix get_Transformed_reference_coordinates()
		{
			return transformed_subset_REF_PDB_coordinates;
		}

	public List<Integer> get_conformation_outliers()
		{

			return conformation_outliers;
		}

	public Matrix get_conf_Z_scores()
		{

			return conf_Z_scores;
		}

	public int getNumber_of_residues()
		{

			return number_of_residues;
		}

	public double get_z_cut()
		{

			return z_cut_off;
		}

	public double get_percent()
		{

			return percent_cut;
		}

}
