package jed;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

import Jama.Matrix;

/**
 * JED class VIZ_Driver: Driver program to visualize Cartesian PCA modes.
 * Constructs sets of 20 PDB files and a Pymol(TM) Script to animate the PCA modes derived from the Cartesian subset.
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
public class VIZ_Driver
{

	static String line, input_path, out_dir, PDB, eigenvectors, evals, modes, maxes, mins, file_name_head, CA = "CA", C = "C", N = "N", O = "O", date, type = "Test";
	static int number_of_jobs, job_number, number_of_frames, number_of_cycles, number_Of_Input_Lines, line_count, number_of_modes_viz, number_of_residues, ROWS_Evectors, ROWS_Modes, COLS;
	static double mode_amplitude, normed_mode_amplitude;
	static final double FLOOR = 1.00E-3, delta_y = 99, pi = Math.PI;
	static List<Double> eigenvalues, pca_mode_maxs, pca_mode_mins;
	static Matrix top_evectors, top_square_pca_modes;
	static BufferedWriter output_file_writer;
	static Vector<Atom> atoms, atoms_essential;
	static File log;
	static BufferedReader input_reader;
	static BufferedWriter log_writer;
	static DateUtils now;
	static StringTokenizer sToken;
	static NumberFormat nf;
	static RoundingMode rm;
	static long startTime, endTime, totalTime;
	static List<String> lines;
	static FortranFormat formatter;
	static PDB_File_Parser parser;
	static boolean do_individual, exist, success, OK;

	private static void read_data_files()
		{
			atoms = PDB_IO.Read_PDB(PDB);
			atoms_essential = PDB_IO.Read_PDB(PDB);
			top_evectors = Matrix_IO.read_Matrix(eigenvectors);
			top_square_pca_modes = Matrix_IO.read_Matrix(modes);
			eigenvalues = List_IO.read_List(evals, "Double");
			pca_mode_maxs = List_IO.read_List(maxes, "Double");
			pca_mode_mins = List_IO.read_List(mins, "Double");
		}

	private static void get_Mode_Visualizations()
		{

			ROWS_Evectors = top_evectors.getRowDimension();
			ROWS_Modes = top_square_pca_modes.getRowDimension();
			number_of_residues = ROWS_Modes;
			COLS = top_evectors.getColumnDimension();
			file_name_head = out_dir + "ss_" + number_of_residues;
			formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
			formatter.setAddReturn(true);
			parser = new PDB_File_Parser();

			for (int outer = 0; outer < number_of_modes_viz; outer++) // iterates over the modes
				{

					double MODE_MIN = pca_mode_mins.get(outer);
					double MODE_MAX = pca_mode_maxs.get(outer);
					if (MODE_MIN < FLOOR) MODE_MIN = FLOOR;
					double LOG_MODE_MIN = Math.log10(MODE_MIN);
					double LOG_MODE_MAX = Math.log10(MODE_MAX);
					double delta_x = (LOG_MODE_MAX - LOG_MODE_MIN);
					double slope = (delta_y / delta_x);
					double y_min = (slope * LOG_MODE_MIN);
					String f_index = "";
					int frame_index = 0;
					Matrix evector = top_evectors.getMatrix(0, ROWS_Evectors - 1, outer, outer);
					Matrix mode = top_square_pca_modes.getMatrix(0, ROWS_Modes - 1, outer, outer);
					double[] evector_array = evector.getColumnPackedCopy();
					Arrays.sort(evector_array);
					double evector_max = Math.max(Math.abs(evector_array[0]), Math.abs(evector_array[ROWS_Evectors - 1]));
					normed_mode_amplitude = Math.abs(mode_amplitude - evector_max);
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
													double shift_x = v_x * w * normed_mode_amplitude;
													double shift_y = v_y * w * normed_mode_amplitude;
													double shift_z = v_z * w * normed_mode_amplitude;
													a.x = x_coord + shift_x;
													a.y = y_coord + shift_y;
													a.z = z_coord + shift_z;
													double bff = (mode.get(index, 0));
													if (bff < FLOOR) bff = FLOOR;
													double log_bff = Math.log10(bff);
													double bf = ((slope * log_bff) - y_min);
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
							System.err.println("IO Exception thrown. Could not write the mode file: " + file_name_head + "_Mode_" + (outer + 1) + "_frame_" + f_index + ".pdb");
							io.printStackTrace();
							System.exit(0);
						}
					write_Pymol_Script(outer);
				}
		}

	private static void get_Essential_Visualization() throws IOException
		{

			ROWS_Evectors = top_evectors.getRowDimension();
			ROWS_Modes = top_square_pca_modes.getRowDimension();
			number_of_residues = ROWS_Modes;
			COLS = top_evectors.getColumnDimension();
			file_name_head = out_dir + "ss_" + number_of_residues;
			formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
			formatter.setAddReturn(true);
			parser = new PDB_File_Parser();

			/* Get the top combined square mode */
			Matrix sum = top_square_pca_modes.getMatrix(0, ROWS_Modes - 1, 0, 0);
			for (int i = 1; i < number_of_modes_viz; i++)
				{
					Matrix plus = top_square_pca_modes.getMatrix(0, ROWS_Modes - 1, i, i);
					sum = sum.plus(plus);
				}
			// sum.print(12, 3);

			/* Norm and sort the vector */
			Matrix sum_normed = Projector.get_Normed_array(sum);
			double[] sorted_sum_normed = sum_normed.getColumnPackedCopy();
			Arrays.sort(sorted_sum_normed);

			// sum_normed.print(12, 3);

			/* Establish the Log coloring scheme */
			double MODE_MIN = sorted_sum_normed[0];
			double MODE_MAX = sorted_sum_normed[COLS - 1];
			if (MODE_MIN < FLOOR) MODE_MIN = FLOOR;
			double LOG_MODE_MIN = Math.log10(MODE_MIN);
			double LOG_MODE_MAX = Math.log10(MODE_MAX);
			double delta_x = (LOG_MODE_MAX - LOG_MODE_MIN);
			double slope = (delta_y / delta_x);
			double y_min = (slope * LOG_MODE_MIN);
			final double eigenvalue_max = eigenvalues.get(0); // the largest eigenvalue sets the wave number for the first mode
			System.out.println("Combined Mode MAX = " + nf.format(MODE_MAX));
			System.out.println("Combined Mode MIN = " + nf.format(MODE_MIN));
			System.out.println("Eigenvalue max = " + nf.format(eigenvalue_max));
			System.out.println("The number of cycles is: " + number_of_cycles);
			System.out.println("The number of frames to generate is: " + number_of_frames);
			System.out.println("The mode amplitude is: " + mode_amplitude);
			int frame_index = 0;
			String f_index = "";

			for (int t = 0; t < number_of_frames; t++) // loop for perturbing the eigenvector components sinusoidally over number of frames: FRAME LOOP
				{
					f_index = String.format("%03d", frame_index + 1);
					frame_index++;
					normed_mode_amplitude = 1.00;

					double omega = 0;
					double weight = 0;

					for (int mode = 0; mode < number_of_modes_viz; mode++) // iterates over the modes: MODE LOOP
						{
							Matrix evector = top_evectors.getMatrix(0, ROWS_Evectors - 1, mode, mode);
							double eval = eigenvalues.get(mode);
							double A_k = Math.sqrt(eval / eigenvalue_max);
							omega = ((2 * number_of_cycles * pi / number_of_frames) * Math.sqrt(eigenvalue_max / eval)); // set the frequency
							weight = A_k * Math.sin(omega * t); // sine function ensures harmonic motion where the first structure is unperturbed;
							int index = 0;
							int count = 0;
							for (Atom a : atoms_essential) // Iterates through the vector of atoms; shifts all backbone atoms along the eigenvector: ATOM LOOP
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
											double bff = (sum_normed.get(index, 0));
											if (bff < FLOOR) bff = FLOOR;
											double log_bff = Math.log10(bff);
											double bf = ((slope * log_bff) - y_min);
											a.setB_factor(bf);
										} else
										{
											a.setB_factor(0);
										}
								}
						}
					String output_file = file_name_head + "_Essential_Modes_" + (number_of_modes_viz) + "_frame_" + f_index + ".pdb";
					output_file_writer = new BufferedWriter(new FileWriter(output_file));
					parser.write_PDB(output_file_writer, atoms_essential, formatter);
					output_file_writer.close();
					write_Pymol_Script_Essential(number_of_modes_viz);
				}
		}

	private static void write_Pymol_Script(int mode_number)
		{
			try
				{
					String name = "ss_" + number_of_residues;
					File pymol_script_file = new File(file_name_head + "_Mode_" + (mode_number + 1) + ".pml");
					BufferedWriter script_file_writer = new BufferedWriter(new FileWriter(pymol_script_file));
					script_file_writer.write("from pymol import cmd" + "\n");
					script_file_writer.write("from pymol.cgo import *" + "\n");
					script_file_writer.write("bg_color white" + "\n");
					script_file_writer.write("from glob import glob" + "\n");
					script_file_writer.write("filelist = glob (" + " \"" + name + "_Mode_" + (mode_number + 1) + "*.pdb\" )" + "\n");
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
					System.err.println("IOException thrown. Could not write the Pymol script file: " + file_name_head + "_Mode_" + (mode_number + 1) + ".pml");
					io.getMessage();
					io.getStackTrace();
				}
		}

	private static void write_Pymol_Script_Essential(int modes)
		{
			try
				{
					String name = "ss_" + number_of_residues;
					File pymol_script_file = new File(file_name_head + "_Essential_Modes_" + (modes) + ".pml");
					BufferedWriter script_file_writer = new BufferedWriter(new FileWriter(pymol_script_file));
					script_file_writer.write("from pymol import cmd" + "\n");
					script_file_writer.write("from pymol.cgo import *" + "\n");
					script_file_writer.write("bg_color white" + "\n");
					script_file_writer.write("from glob import glob" + "\n");
					script_file_writer.write("filelist = glob (" + " \"" + name + "_Essential_Modes_" + (modes) + "*.pdb\" )" + "\n");
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
					System.err.println("IOException thrown. Could not write the Pymol script file: " + file_name_head + "_Mode_" + (modes) + ".pml");
					io.getMessage();
					io.getStackTrace();
				}
		}

	private static void read_input_file()
		{
			number_Of_Input_Lines = 0;
			line_count = 0;
			lines = new ArrayList<String>();
			System.out.println("Below is the input file that was read: " + input_path);
			System.out.println("--------------------------------------------------------------------------------------------------------------------------------");
			try
				{
					while ((line = input_reader.readLine()) != null && line.length() >= 1)
						{
							lines.add(line);
							System.out.println(line);
							number_Of_Input_Lines++;
						}

					input_reader.close();
					System.out.println("--------------------------------------------------------------------------------------------------------------------------------");
					System.out.println("The number of lines of parameters in the input file is: " + number_Of_Input_Lines + "\n");
					if (number_Of_Input_Lines < 3)
						{
							System.err.println("INSUFFICIENT DATA IN THE INPUT FILE:");
							System.err.println("THERE MUST BE AT LEAST 3 LINES OF PARAMETERS.");
							System.err.println("Terminating program execution.");
							System.exit(0);
						}

				} catch (IOException e)
				{
					System.err.println("IOException thrown. Could not read the input file. Program will terminate.\n");
					e.printStackTrace();
					System.exit(0);
				}
			System.gc();
		}

	private static void read_batch_parameters()
		{
			line_count = 0;
			System.out.println("Reading line " + (line_count + 1)); // Reads line 1, the number of jobs
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			String test = sToken.nextToken();
			OK = Test_Numeric_Type.test_Integer(test);
			if (!OK || Integer.parseInt(test) < 1)
				{
					System.err.println("Expected Number of Jobs to be a positive integer, but got: " + test);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			if (OK) number_of_jobs = Integer.parseInt(test);
			System.out.println("\tThe number of jobs =  " + number_of_jobs);
			line_count++;
			/* *********************************************************************************************************************************** */
			System.out.println("Reading line " + (line_count + 1)); // Reads the divider line between the number of jobs and the first job
			line = lines.get(line_count);
			System.out.println(line + "\n");
			line_count++;
			/* *********************************************************************************************************************************** */

		}

	private static void read_job_parameters()
		{
			if (line_count >= lines.size())
				{
					System.err.println("No more data in input file for remaining jobs in batch.");
					System.err.println("User specified too many jobs.");
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			/* ************************************************************************************************************************************ */
			System.out.println("Reading line " + (line_count + 1) + ", the FIRST LINE of job: " + (job_number + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			String test = sToken.nextToken();
			OK = Test_Numeric_Type.test_Integer(test);
			if (OK) number_of_modes_viz = Integer.parseInt(test);
			System.out.println("\tThe number of modes to visualize = " + number_of_modes_viz);
			if (!OK || number_of_modes_viz < 1)
				{
					System.err.println("Expected the number of modes to visualize to be a positive integer, but got: " + test);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			test = sToken.nextToken();
			OK = Test_Numeric_Type.test_Double(test);
			if (OK) mode_amplitude = Double.parseDouble(test);
			System.out.println("\tThe Mode Amplitude = " + mode_amplitude);
			if (!OK || mode_amplitude < 0)
				{
					System.err.println("Expected Mode Amplitude to be a positive decimal, but got: " + test);
					System.err.println("Setting Mode Amplitude to default value of 1.5");
					mode_amplitude = 1.5;
				}
			test = sToken.nextToken();
			OK = Test_Numeric_Type.test_Integer(test);
			if (OK) number_of_frames = Integer.parseInt(test);
			System.out.println("\tThe Number of frames to generate = " + number_of_frames);
			if (!OK || number_of_frames < 0)
				{
					System.err.println("Expected Number of frames to be a positive integer, but got: " + test);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			test = sToken.nextToken();
			OK = Test_Numeric_Type.test_Integer(test);
			if (OK) number_of_cycles = Integer.parseInt(test);
			System.out.println("\tThe Number of Cycles to generate = " + number_of_cycles);
			if (!OK || number_of_cycles < 0)
				{
					System.err.println("Expected Number of Cycles to be a positive integer, but got: " + test);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			test = sToken.nextToken();
			if (test == "1") do_individual = true;
			line_count++;
			/* ************************************************************************************************************************************ */
			System.out.println("Reading line " + (line_count + 1) + ", the SECOND LINE of job: " + (job_number + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			PDB = sToken.nextToken();
			System.out.println("\tPath to the Reference PDB file = " + PDB);
			exist = new File(PDB).exists();
			if (!exist)
				{
					System.err.println("The file does not exist:  " + PDB);
					System.exit(0);
				}
			line_count++;
			/* ************************************************************************************************************************************ */
			System.out.println("Reading line " + (line_count + 1) + ", the THIRD LINE of job: " + (job_number + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			evals = sToken.nextToken();
			System.out.println("\tPath to the Eigenvalues file = " + evals);
			exist = new File(evals).exists();
			if (!exist)
				{
					System.err.println("The file does not exist:  " + evals);
					System.exit(0);
				}
			line_count++;
			/* ************************************************************************************************************************************ */
			System.out.println("Reading line " + (line_count + 1) + ", the FOURTH LINE of job: " + (job_number + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			eigenvectors = sToken.nextToken();
			System.out.println("\tPath to the Eigenvector file = " + eigenvectors);
			exist = new File(eigenvectors).exists();
			if (!exist)
				{
					System.err.println("The file does not exist:  " + eigenvectors);
					System.exit(0);
				}
			line_count++;
			/* ************************************************************************************************************************************ */
			System.out.println("Reading line " + (line_count + 1) + ", the FIFTH LINE of job: " + (job_number + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			modes = sToken.nextToken();
			System.out.println("\tPath to the Modes file = " + modes);
			exist = new File(modes).exists();
			if (!exist)
				{
					System.err.println("The file does not exist:  " + modes);
					System.exit(0);
				}
			line_count++;
			/* ************************************************************************************************************************************ */
			System.out.println("Reading line " + (line_count + 1) + ", the SIXTH LINE of job: " + (job_number + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			maxes = sToken.nextToken();
			System.out.println("\tPath to the Mode MAXES file = " + maxes);
			exist = new File(maxes).exists();
			if (!exist)
				{
					System.err.println("The file does not exist:  " + maxes);
					System.exit(0);
				}
			line_count++;
			/* ************************************************************************************************************************************ */
			System.out.println("Reading line " + (line_count + 1) + ", the SEVENTH LINE of job: " + (job_number + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			mins = sToken.nextToken();
			System.out.println("\tPath to the Mode MINS file = " + mins);
			exist = new File(mins).exists();
			if (!exist)
				{
					System.err.println("The file does not exist:  " + mins);
					System.exit(0);
				}
			line_count++;
			/* ************************************************************************************************************************************ */
			System.out.println("Reading line " + (line_count + 1) + ", the EIGHTH LINE of job: " + (job_number + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			out_dir = sToken.nextToken();
			System.out.println("\tOutput Directory = " + out_dir);
			if (!(out_dir.endsWith("/") || out_dir.endsWith("\\")))
				{
					System.err.println("Expected the Output Directory to end with / or \\\\, but got: " + line);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			exist = new File(out_dir).exists();
			if (!exist)
				{
					System.err.println("\tThe output directory does not exist.");
					System.err.println("\tAttempting to create it:");
					boolean success = (new File(out_dir)).mkdirs();
					if (success) System.err.println("\t\tSuccess.");
					if (!success)
						{
							System.err.println("Failed to create the output directory:  " + out_dir);
							System.exit(0);
						}
				}
			line_count++;
			/* ************************************************************************************************************************************ */
			line = lines.get(line_count); // Reads the divider line between jobs
			System.out.println("Reading line " + (line_count + 1));
			System.out.println(line + "\n");
			line_count++;
			/* ************************************************************************************************************************************ */
		}

	private static void initialize_Job_Log()
		{
			try
				{
					log = new File(out_dir + "VIZ_LOG.txt");
					log_writer = new BufferedWriter(new FileWriter(log));
					log_writer.write("JED Cartesian Mode Visualization Driver version 1.0" + "\n");
					log_writer.write("Reference PDB file: " + PDB + "\n");
					log_writer.write("Eigenvalues file: " + evals + "\n");
					log_writer.write("Eigenvectors file: " + eigenvectors + "\n");
					log_writer.write("Modes file: " + modes + "\n");
					log_writer.write("Mode Maxes file: " + maxes + "\n");
					log_writer.write("Mode Mins file: " + mins + "\n");
					log_writer.write("Output directory: " + out_dir + "\n");
					log_writer.write("Number of Frames = " + (number_of_frames) + "\n");
					log_writer.write("Number of Cycles = " + (number_of_cycles) + "\n");
					log_writer.write("Mode amplitutde = " + nf.format(mode_amplitude) + "\n");

					log_writer.flush();

				} catch (IOException e)
				{
					System.err.println("Could not write the VIZ_LOG file: " + out_dir + "VIZ_LOG.txt");
					System.err.println("Program terminating.\n");
					e.printStackTrace();
					System.exit(0);
				}
			System.gc();
		}

	private static void write_VIZ_Log()
		{
			try
				{
					log_writer.write("\nPerforming Cartesian Mode Visualization on Top  " + number_of_modes_viz + "  individual cPCA modes.\n");
					log_writer.write("MODE AMPLITUDE = " + nf.format(mode_amplitude) + "\n");
					log_writer.flush();

				} catch (IOException e)
				{
					System.err.println("Could not write the VIZ_LOG file: " + out_dir + "VIZ_LOG.txt");
					System.err.println("Program terminating.\n");
					e.printStackTrace();
					System.exit(0);
				}
			System.gc();
		}

	private static void write_Essential_VIZ_Log()
		{
			try
				{
					log_writer.write("\nPerforming Cartesian Essential Mode Visualization, combining the top  " + number_of_modes_viz + "  cPCA modes.\n");
					log_writer.flush();

				} catch (IOException e)
				{
					System.err.println("Could not write the VIZ_LOG file: " + out_dir + "VIZ_LOG.txt");
					System.err.println("Program terminating.\n");
					e.printStackTrace();
					System.exit(0);
				}
			System.gc();
		}

	public static void main(String[] args) throws IOException
		{

			rm = RoundingMode.HALF_UP;
			nf = NumberFormat.getInstance();
			nf.setMaximumFractionDigits(3);
			nf.setMinimumFractionDigits(3);
			nf.setRoundingMode(rm);

			System.out.println("Running VIZ Driver: ");
			System.out.println("Getting the input file: ");
			try
				{

					input_path = "viz.txt";
					String WD = System.getProperty("user.dir");
					String in_path = WD + File.separator + input_path;

					if (args.length >= 1)
						{
							input_path = args[0];
							System.out.println("The path to the input file must be the first argument:");
							System.out.println("These are the command args:");
							for (int i = 0; i < args.length; i++)
								{
									System.out.println("Arg " + (i + 1) + " Value = " + args[i]);
								}
							in_path = input_path;
						}
					System.out.println("Working Directory = " + WD);
					System.out.println("Input File Path = " + in_path);
					boolean check = new File(in_path).exists();
					if (!check)
						{
							System.err.println("The entered Input File does not exist: " + in_path);
							System.err.println("Terminating program execution.");
							System.exit(0);
						}
					input_reader = new BufferedReader(new FileReader(in_path));

				} catch (FileNotFoundException e)
				{
					System.err.println("Could not find the input file: " + input_path);
					System.err.println("Program terminating.\n");
					e.printStackTrace();
					System.exit(0);
				}

			System.out.println("Reading Input File... ");

			read_input_file();
			read_batch_parameters();

			for (job_number = 0; job_number < number_of_jobs; job_number++)
				{
					System.out.println("Now running job " + (job_number + 1) + " of " + number_of_jobs);

					read_job_parameters();
					initialize_Job_Log();
					read_data_files();

					if (do_individual)
						{
							System.out.println("Performing Cartesian Mode Visualizations... ");
							startTime = System.nanoTime();
							get_Mode_Visualizations();
							write_VIZ_Log();
							endTime = System.nanoTime();
						}

					System.out.println("Performing Essential Cartesian Mode Visualization... ");
					startTime = System.nanoTime();
					get_Essential_Visualization();
					write_Essential_VIZ_Log();
					endTime = System.nanoTime();

					totalTime = endTime - startTime;

					System.out.println("Done. (" + nf.format(totalTime / 1000000000.0) + " seconds)");

					date = DateUtils.now();

					try
						{
							log_writer.write("\nCartesian Mode Visualization Completed: " + date);
							log_writer.close();

						} catch (IOException e)
						{
							System.err.println("Could not write the VIZ_LOG file");
							e.printStackTrace();
						}

					System.out.println("JED Cartesian Mode Visualization Driver successfully completed: " + date);
				}
		}
}
