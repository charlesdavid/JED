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
	final double FLOOR = 1.00E-3, delta_y = 99, pi = Math.PI;
	double[] pca_mode_maxs, pca_mode_mins;
	List<Double> eigenvalues;
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
	 * @param residues
	 *            The list of original residues from the PDB file
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
	JED_PCA_Mode_Vizualization(String dir, String des, Vector<Atom> atms, List<Double> evals, Matrix evects, Matrix modes, double[] maxs, double[] mins, double mode_amp, int num_modes, String T)
		{
			directory = dir;
			description = des;
			type = T;
			out_dir = directory + "JED_RESULTS_" + description + "/" + "VIZ/" + type + "/";
			exist = new File(out_dir).exists();
			if (!exist) success = (new File(out_dir)).mkdirs();
			atoms = atms;
			eigenvalues = evals;
			top_evectors = evects;
			square_pca_modes = modes;
			ROWS_Evectors = evects.getRowDimension();
			ROWS_Modes = square_pca_modes.getRowDimension();
			number_of_residues = ROWS_Modes;
			COLS = evects.getColumnDimension();
			number_of_modes_viz = num_modes;
			pca_mode_maxs = maxs;
			pca_mode_mins = mins;
			mode_amplitude = mode_amp;
			formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
			formatter.setAddReturn(true);
			parser = new PDB_File_Parser();
			file_name_head = out_dir + "ss_" + number_of_residues;
		}

	/**
	 * Method that generates the PDB files and the PyMOL (TM) scripts for the mode visualization
	 */
	void get_Mode_Visualizations_SS()
		{

			ROWS_Evectors = top_evectors.getRowDimension();
			ROWS_Modes = square_pca_modes.getRowDimension();
			number_of_residues = ROWS_Modes;
			COLS = top_evectors.getColumnDimension();
			file_name_head = out_dir + "ss_" + number_of_residues;
			formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
			formatter.setAddReturn(true);
			parser = new PDB_File_Parser();

			for (int outer = 0; outer < number_of_modes_viz; outer++) // iterates over the modes
				{

					double MODE_MIN = pca_mode_mins[outer];
					double MODE_MAX = pca_mode_maxs[outer];
					if (MODE_MIN < FLOOR) MODE_MIN = FLOOR;
					double LOG_MODE_MIN = Math.log10(MODE_MIN);
					double LOG_MODE_MAX = Math.log10(MODE_MAX);
					double delta_x = (LOG_MODE_MAX - LOG_MODE_MIN);
					double slope = (delta_y / delta_x);
					String f_index = "";
					int frame_index = 0;
					Matrix evector = top_evectors.getMatrix(0, ROWS_Evectors - 1, outer, outer);
					Matrix mode = square_pca_modes.getMatrix(0, ROWS_Modes - 1, outer, outer);
					try
						{
							for (double mc = 0; mc < (2 * Math.PI); mc += (Math.PI / 10)) // loop for perturbing the eigenvector components sinusoidally
								{
									f_index = String.format("%03d", frame_index + 1);
									String output_file = file_name_head + "_Mode_" + (outer + 1) + "_" + type + "_frame_" + f_index + ".pdb";
									output_file_writer = new BufferedWriter(new FileWriter(output_file));
									int index = 0;
									int count = 0;
									for (Atom a : atoms) // Iterates through the vector of atoms; shifts all backbone atoms along the eigenvector
										{
											index = (count / 4);  // there are 4 atom symbols in the if statement
											if (a.symbol.equals(N) || a.symbol.equals(CA) || a.symbol.equals(C) || a.symbol.equals(O))
												{
													count++;
													double x_coord = a.x;
													double y_coord = a.y;
													double z_coord = a.z;
													double v_x = evector.get(index, 0);
													double v_y = evector.get((index + number_of_residues), 0);
													double v_z = evector.get((index + 2 * number_of_residues), 0);
													double w = Math.sin(mc); // sine function ensures the first structure is unperturbed;
																			 // preserves chain connectivity in PyMol movies...
													double shift_x = v_x * w * mode_amplitude;
													double shift_y = v_y * w * mode_amplitude;
													double shift_z = v_z * w * mode_amplitude;
													a.x = x_coord + shift_x;
													a.y = y_coord + shift_y;
													a.z = z_coord + shift_z;
													double bff = (mode.get(index, 0));
													if (bff < MODE_MIN) bff = MODE_MIN;
													if (bff > MODE_MAX) bff = MODE_MAX;
													double BFF = bff / MODE_MIN;
													double log_bff = Math.log10(BFF);
													double bf = (slope * log_bff);
													a.setB_factor(bf);
												} else
												{
													a.setB_factor(0);
												}
										}
									parser.write_PDB(output_file_writer, atoms, formatter);
									output_file_writer.close();
									frame_index++;
									System.gc();

								}

						} catch (IOException io)
						{
							System.err.println("IO Exception thrown. Could not write the mode file: " + output_file);
							io.printStackTrace();
							System.exit(0);
						}
					write_Pymol_Script(outer);
				}
		}

	void get_Essential_Visualization()
		{

			ROWS_Evectors = top_evectors.getRowDimension();
			ROWS_Modes = square_pca_modes.getRowDimension();
			number_of_residues = ROWS_Modes;
			COLS = top_evectors.getColumnDimension();
			file_name_head = out_dir + "ss_" + number_of_residues;
			formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
			formatter.setAddReturn(true);
			parser = new PDB_File_Parser();

			/* Get the top combined square mode */
			Matrix sum = square_pca_modes.getMatrix(0, ROWS_Modes - 1, 0, 0);
			double eval_abs = Math.abs(eigenvalues.get(0));
			Matrix w_sum = sum.times(eval_abs);
			for (int i = 1; i < number_of_modes_viz; i++)
				{
					Matrix plus = square_pca_modes.getMatrix(0, ROWS_Modes - 1, i, i);
					eval_abs = Math.abs(eigenvalues.get(i));
					Matrix w_plus = plus.times(eval_abs);
					w_sum = w_sum.plus(w_plus);
				}

			/* Norm and sort the vector */
			Matrix w_sum_normed = Projector.get_Normed_array(w_sum);
			double[] sorted_w_sum_normed = w_sum_normed.getColumnPackedCopy();
			Arrays.sort(sorted_w_sum_normed);

			/* Establish the Log coloring scheme */
			double MODE_MIN = sorted_w_sum_normed[0];
			double MODE_MAX = sorted_w_sum_normed[ROWS_Modes - 1];
			if (MODE_MIN < FLOOR) MODE_MIN = FLOOR;
			double LOG_MODE_MIN = Math.log10(MODE_MIN);
			double LOG_MODE_MAX = Math.log10(MODE_MAX);
			double delta_x = (LOG_MODE_MAX - LOG_MODE_MIN);
			double slope = (delta_y / delta_x);

			final double eigenvalue_max = eigenvalues.get(0); // the largest eigenvalue sets the wave number for the first mode
			int frame_index = 0;
			String f_index = "";
			int number_of_frames = 20; // hard code for 20 frames
			int number_of_cycles = 2; // hard code for 2 cycles of the fundamental mode
			int num_of_essential_modes = 5; // hard code for using only the top 5 modes
			if (number_of_modes_viz < num_of_essential_modes) num_of_essential_modes = number_of_modes_viz;
			for (int t = 0; t < number_of_frames; t++) // loop for perturbing the eigenvector components sinusoidally over number of frames: FRAME LOOP
				{
					f_index = String.format("%03d", frame_index + 1);
					frame_index++;
					double omega = 0;
					double weight = 0;

					for (int mode = 0; mode < num_of_essential_modes; mode++) // iterates over the modes: MODE LOOP
						{
							Matrix evector = top_evectors.getMatrix(0, ROWS_Evectors - 1, mode, mode);
							double eval = eigenvalues.get(mode);
							double A_k = Math.sqrt(eval / eigenvalue_max);
							omega = ((2 * number_of_cycles * pi / number_of_frames) * Math.sqrt(eigenvalue_max / eval));
							weight = A_k * Math.sin(omega * t); // sine function ensures harmonic motion where the first structure is unperturbed;
							if (type.equals("PCORR")) // adjust for inverse
								{
									omega = ((2 * number_of_cycles * pi / number_of_frames) * Math.sqrt(eval / eigenvalue_max)); // set the frequency for PCORR
									A_k = Math.sqrt(eigenvalue_max / eval);
									weight = A_k * Math.sin(omega * t);
								}
							int index = 0;
							int count = 0;
							for (Atom a : atoms) // Iterates through the vector of atoms; shifts all backbone atoms along the eigenvector: ATOM LOOP
								{
									index = (count / 4);  // there are 4 atom symbols in the if statement
									if (a.symbol.equals(N) || a.symbol.equals(CA) || a.symbol.equals(C) || a.symbol.equals(O))
										{
											count++;
											double x_coord = a.x;
											double y_coord = a.y;
											double z_coord = a.z;
											double v_x = evector.get(index, 0);
											double v_y = evector.get((index + number_of_residues), 0);
											double v_z = evector.get((index + 2 * number_of_residues), 0);
											double shift_x = v_x * weight * mode_amplitude;
											double shift_y = v_y * weight * mode_amplitude;
											double shift_z = v_z * weight * mode_amplitude;
											a.x = x_coord + shift_x;
											a.y = y_coord + shift_y;
											a.z = z_coord + shift_z;
											double bff = (w_sum_normed.get(index, 0));
											if (bff < MODE_MIN) bff = MODE_MIN;
											if (bff > MODE_MAX) bff = MODE_MAX;
											double BFF = (bff / MODE_MIN);
											double log_bff = Math.log10(BFF);
											double bf = (slope * log_bff);
											a.setB_factor(bf);

										} else
										{
											a.setB_factor(0);
										}
								}
						}
					String output_file = file_name_head + "_Essential_Modes_" + number_of_modes_viz + "_" + type + "_frame_" + f_index + ".pdb";

					try
						{
							output_file_writer = new BufferedWriter(new FileWriter(output_file));
							parser.write_PDB(output_file_writer, atoms, formatter);
							output_file_writer.close();
							write_Pymol_Script_Essential(number_of_modes_viz);
						} catch (IOException io)
						{
							System.err.println("IO Exception thrown. Could not write the mode file: " + output_file);
							io.printStackTrace();
							System.exit(0);
						}

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
					script_file_writer.write("filelist = glob (" + " \"" + name + "_Mode_" + (mode_number + 1) + "_" + type + "_frame*.pdb\" )" + "\n");
					script_file_writer.write("for file in filelist: cmd.load( file, " + "\"" + "Mode_" + (mode_number + 1) + "\" )" + "\n");
					script_file_writer.write("hide lines, " + "Mode_" + (mode_number + 1) + "\n");
					script_file_writer.write("show cartoon, " + "Mode_" + (mode_number + 1) + "\n");
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
					System.err.println("IOException thrown. Could not write the PyMOL (TM) script file: " + file_name_head + "_Mode_" + (mode_number + 1) + "_" + type + ".pml");
					io.getMessage();
					io.getStackTrace();
				}
		}

	private void write_Pymol_Script_Essential(int modes)
		{
			String name = "ss_" + number_of_residues;
			String filename = file_name_head + "_Essential_Modes_" + number_of_modes_viz + "_" + type + ".pml";
			try
				{
					File pymol_script_file = new File(filename);
					BufferedWriter script_file_writer = new BufferedWriter(new FileWriter(pymol_script_file));
					script_file_writer.write("from pymol import cmd" + "\n");
					script_file_writer.write("from pymol.cgo import *" + "\n");
					script_file_writer.write("bg_color white" + "\n");
					script_file_writer.write("from glob import glob" + "\n");
					script_file_writer.write("filelist = glob (" + " \"" + name + "_Essential_Modes_" + number_of_modes_viz + "_" + type + "_frame*.pdb\" )" + "\n");
					script_file_writer.write("for file in filelist: cmd.load( file, " + "\"" + "Essential_Modes_" + (modes) + "\" )" + "\n");
					script_file_writer.write("hide lines, " + "Essential_Modes_" + (modes) + "\n");
					script_file_writer.write("show cartoon, " + "Essential_Modes_" + (modes) + "\n");
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
					System.err.println("IOException thrown. Could not write the PyMOL(TM) script file: " + filename);
					io.getMessage();
					io.getStackTrace();
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
