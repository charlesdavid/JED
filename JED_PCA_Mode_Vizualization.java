package jed;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Vector;

import Jama.Matrix;

/**
 * JED class JED_PCA_Mode_Visualization: Constructs sets of 20 PDB files and a Pymol(TM) Script to animate the PCA modes derived from the Cartesian subset.
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
 * along with this program. If not, see <http://www.gnu.org/license>.
 *
 * @author Dr.Charles David
 */

/**
 * @author Charles
 *
 */
public class JED_PCA_Mode_Vizualization
	{

		String directory, out_dir, description, output_file, file_name_head, type, CA = "CA", C = "C", N = "N", O = "O";
		int number_of_modes_viz, number_of_residues, ROWS_Evectors, ROWS_Modes, COLS;
		double mode_amplitude;
		final double FLOOR = 1.00E-3, delta_y = 99;
		double[] pca_mode_max, pca_mode_min;
		List<Integer> residue_list_original;
		Matrix top_evectors, square_pca_modes;
		BufferedWriter output_file_writer;
		Vector<Atom> atoms;
		FortranFormat formatter;
		PDB_File_Parser parser;
		boolean exist, success;

		/**
		 *
		 * Constructor for generating the sets of PDB files for visualizing the Cartesian PCA modes.
		 *
		 * @param dir
		 *            The working directory
		 * @param des
		 *            The job description
		 * @param atms
		 *            The atoms that comprise the subset
		 * @param evects
		 *            The eigenvectors for the top cPCA modes
		 * @param modes
		 *            The top cPCA modes
		 * @param maxs
		 *            The array containing the maximum component of each cPCA mode
		 * @param mins
		 *            The array containing the minimum component of each cPCA mode
		 * @param mode_amp
		 *            The mode amplitude, which controls how far the CA, C, N, and O atoms are displaced
		 * @param num_modes
		 *            The number of cPCA modes for which to generate structures
		 * @param T
		 *            The PCA model: COV, CORR, or PCORR
		 */
		JED_PCA_Mode_Vizualization(String dir, String des, List<Integer> residues, Vector<Atom> atms, Matrix evects, Matrix modes, double[] maxs, double[] mins,
				double mode_amp, int num_modes, String T)
			{

				directory = dir;
				description = des;
				residue_list_original = residues;
				type = T;
				out_dir = directory + "JED_RESULTS_" + description + "/" + "VIZ/" + type + "/";
				exist = new File(out_dir).exists();
				if (!exist) success = (new File(out_dir)).mkdirs();
				atoms = atms;
				top_evectors = evects;
				square_pca_modes = modes;
				ROWS_Evectors = evects.getRowDimension();
				ROWS_Modes = square_pca_modes.getRowDimension();
				number_of_residues = ROWS_Modes;
				COLS = evects.getColumnDimension();
				number_of_modes_viz = num_modes;

				if (number_of_modes_viz > COLS)
					{
						System.err.println("FATAL ERROR!");
						System.err.println("Number of Cartesian Modes to Visualize REQUESTED: " + number_of_modes_viz);
						System.err.println("Number of Cartesian Modes AVAILABLE: " + COLS);
						System.err.println("The Possible number of Cartesial Modes to Visualize is always <= Number of Cartesian Modes requested.");
						System.err.println("Terminating program execution.");
						System.exit(0);
					}

				pca_mode_max = maxs;
				pca_mode_min = mins;
				mode_amplitude = mode_amp;
				formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
				formatter.setAddReturn(true);
				parser = new PDB_File_Parser();
				file_name_head = out_dir + "ss_" + number_of_residues;

			}

		/**
		 * Method that generates the PDB files and the Pymol (TM) scripts for the mode visualization
		 */
		void get_Mode_Visualizations_SS()
			{

				for (int outer = 0; outer < number_of_modes_viz; outer++) // iterates over the modes
					{

						double MODE_MIN = pca_mode_min[outer];
						double MODE_MAX = pca_mode_max[outer];

						if (MODE_MIN < FLOOR) MODE_MIN = FLOOR;

						double LOG_MODE_MIN = Math.log10(MODE_MIN);
						double LOG_MODE_MAX = Math.log10(MODE_MAX);

						double delta_x = (LOG_MODE_MAX - LOG_MODE_MIN);
						double slope = (delta_y / delta_x);
						double y_min = (slope * LOG_MODE_MIN);

						String f_index = "";
						int frame_index = 0;

						Matrix evector = top_evectors.getMatrix(0, ROWS_Evectors - 1, outer, outer);
						Matrix mode = square_pca_modes.getMatrix(0, ROWS_Modes - 1, outer, outer);

						double[] evector_array = evector.getColumnPackedCopy();
						Arrays.sort(evector_array);
						double evector_max = Math.max(Math.abs(evector_array[0]), Math.abs(evector_array[ROWS_Evectors - 1]));
						mode_amplitude = 1.00 - evector_max;

						try
							{
								for (float mc = -1; mc < 1; mc += .1) // perturbs the eigenvector components
									{
										frame_index = (int) (mc * 10 + 10);

										f_index = String.format("%03d", frame_index + 1);

										String output_file = file_name_head + "_Mode_" + (outer + 1) + "_" + type + "_frame_" + f_index + ".pdb";
										output_file_writer = new BufferedWriter(new FileWriter(output_file));
										for (int index = 0; index < residue_list_original.size(); index++) // iterates over each residue in the subset
											{
												int residue_number = residue_list_original.get(index);
												for (Atom a : atoms) // Iterates through the vector of atoms; shifts all atoms along the eigenvector
													{
														if (a.res_number == residue_number)
															{
																double x_coord = a.x;
																double y_coord = a.y;
																double z_coord = a.z;
																double v_x = evector.get(index, 0);
																double v_y = evector.get((index + number_of_residues), 0);
																double v_z = evector.get((index + 2 * number_of_residues), 0);

																double w_x = v_x * mc * mode_amplitude;
																double w_y = v_y * mc * mode_amplitude;
																double w_z = v_z * mc * mode_amplitude;

																a.x = x_coord + w_x;
																a.y = y_coord + w_y;
																a.z = z_coord + w_z;

																double bff = (mode.get(index, 0));
																if (bff < FLOOR) bff = FLOOR;
																double log_bff = Math.log10(bff);
																double bf = ((slope * log_bff) - y_min);
																a.setB_factor(bf);
															}
													}
											}
										parser.write_PDB(output_file_writer, atoms, formatter);
										output_file_writer.close();
										output_file_writer = null;
										System.gc();
									}

							} catch (IOException io)
							{
								System.err.println("IO Exception thrown. Could not write the mode file: " + file_name_head + "_Mode_" + (outer + 1) + "_frame_"
										+ f_index + ".pdb");
								io.printStackTrace();
								System.exit(0);
							}
						write_Pymol_Script(outer);
					}
			}

		private void write_Pymol_Script(int mode_number)
			{
				try
					{
						String name = "ss_" + number_of_residues;
						File pymol_script_file = new File(file_name_head + "_Mode_" + (mode_number + 1) + "_" + type + ".pml");
						BufferedWriter script_file_writer = new BufferedWriter(new FileWriter(pymol_script_file));
						script_file_writer.write("from pymol import cmd" + "\n");
						script_file_writer.write("from pymol.cgo import *" + "\n");
						script_file_writer.write("bg_color white" + "\n");
						script_file_writer.write("from glob import glob" + "\n");
						script_file_writer
								.write("filelist = glob (" + " \"" + name + "_Mode_" + (mode_number + 1) + "_" + type + "_frame_" + "*.pdb\" )" + "\n");
						script_file_writer.write("for file in filelist: cmd.load( file, " + "\"" + type + "_Mode_" + (mode_number + 1) + "\" )" + "\n");
						script_file_writer.write("hide lines, " + type + "_Mode_" + (mode_number + 1) + "\n");
						script_file_writer.write("show cartoon, " + type + "_Mode_" + (mode_number + 1) + "\n");
						script_file_writer.write("cmd.spectrum(\"b\",selection=\"((all)&*/ca)\",quiet=0)" + "\n");
						script_file_writer.write("orient" + "\n");
						script_file_writer.write("set two_sided_lighting,1" + "\n");
						script_file_writer.write("set cartoon_fancy_helices,1" + "\n");
						script_file_writer.write("set cartoon_highlight_color,grey50" + "\n");
						script_file_writer.write("util.performance(0)" + "\n");
						script_file_writer.write("rebuild" + "\n");
						script_file_writer.close();

					} catch (IOException io)
					{
						System.err.println(
								"IOException thrown. Could not write the Pymol script file: " + file_name_head + "_Mode_" + (mode_number + 1) + ".pml");
						io.printStackTrace();
					}
			}

		double get_Mode_Amplitude()
			{

				return mode_amplitude;
			}

		void set_Mode_Amplitude(double amp)
			{

				mode_amplitude = amp;
			}
	}
