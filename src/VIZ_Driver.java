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

	static String line, input_path, out_dir, PDB, eigenvectors, modes, maxes, mins, file_name_head, CA = "CA", C = "C", N = "N", O = "O", date;
	static int number_Of_Input_Lines, line_count, number_of_modes_viz, number_of_residues, ROWS_Evectors, ROWS_Modes, COLS;
	static double mode_amplitude;
	static final double FLOOR = 1.00E-3, delta_y = 99;
	static int number_of_jobs, job_number;
	static List<Double> pca_mode_maxs, pca_mode_mins;
	static Matrix top_evectors, square_pca_modes;
	static BufferedWriter output_file_writer;
	static Vector<Atom> atoms;
	static File log;
	static BufferedReader input_reader;
	static BufferedWriter log_writer;
	static DateUtils now;
	static StringTokenizer sToken;
	static NumberFormat nf;
	static RoundingMode rm;
	static long startTime, endTime, totalTime;
	static List<String> lines;
	static boolean exist, success, OK;

	private static void read_data_files()
		{
			atoms = PDB_IO.Read_PDB(PDB);
			top_evectors = Matrix_IO.read_Matrix(eigenvectors);
			square_pca_modes = Matrix_IO.read_Matrix(modes);
			pca_mode_maxs = List_IO.read_List(maxes, "Double");
			pca_mode_mins = List_IO.read_List(mins, "Double");
		}

	private static void get_Mode_Visualizations()
		{

			ROWS_Evectors = top_evectors.getRowDimension();
			ROWS_Modes = square_pca_modes.getRowDimension();
			number_of_residues = ROWS_Modes;
			COLS = top_evectors.getColumnDimension();
			file_name_head = out_dir + "ss_" + number_of_residues;

			if (number_of_modes_viz > COLS)
				{
					System.err.println("FATAL ERROR!");
					System.err.println("Number of Cartesian Modes to Visualize REQUESTED: " + number_of_modes_viz);
					System.err.println("Number of Cartesian Modes AVAILABLE: " + COLS);
					System.err.println("The Possible number of Cartesial Modes to Visualize is always <= Number of Cartesian Modes requested.");
					System.err.println("Terminating program execution.");
					System.exit(0);
				}

			for (int outer = 0; outer < number_of_modes_viz; outer++)
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

					try
						{
							for (float mc = -1; mc < 1; mc += .1)
								{
									frame_index = (int) (mc * 10 + 10);
									Matrix evector = top_evectors.getMatrix(0, ROWS_Evectors - 1, outer, outer);
									Matrix d_mode = square_pca_modes.getMatrix(0, ROWS_Modes - 1, outer, outer);

									f_index = String.format("%03d", frame_index + 1);

									String output_file = file_name_head + "_Mode_" + (outer + 1) + "_frame_" + f_index + ".pdb";
									output_file_writer = new BufferedWriter(new FileWriter(output_file));
									int index = 0;
									int count = 0;
									for (Atom a : atoms)
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

													double w_x = v_x * mc * mode_amplitude;
													double w_y = v_y * mc * mode_amplitude;
													double w_z = v_z * mc * mode_amplitude;

													a.x = x_coord + w_x;
													a.y = y_coord + w_y;
													a.z = z_coord + w_z;

													double bff = (d_mode.get(index, 0));
													if (bff < FLOOR) bff = FLOOR;
													double log_bff = Math.log10(bff);

													double bf = ((slope * log_bff) - y_min);

													a.setB_factor(bf);
												} else
												{
													a.setB_factor(0);
												}
										}

									PDB_IO.Write_PDB(output_file, atoms);

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
			/* ************************************************************************************************************************************************** */
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
			/* ************************************************************************************************************************************************** */
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
					System.err.println("Expected Percent to be a positive decimal, but got: " + test);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			line_count++;
			/* ************************************************************************************************************************************************** */
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
			/* ************************************************************************************************************************************************** */
			System.out.println("Reading line " + (line_count + 1) + ", the THIRD LINE of job: " + (job_number + 1));
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
			/* ************************************************************************************************************************************************** */
			System.out.println("Reading line " + (line_count + 1) + ", the FOURTH LINE of job: " + (job_number + 1));
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
			/* ************************************************************************************************************************************************** */
			System.out.println("Reading line " + (line_count + 1) + ", the FIFTH LINE of job: " + (job_number + 1));
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
			/* ************************************************************************************************************************************************** */
			System.out.println("Reading line " + (line_count + 1) + ", the SIXTH LINE of job: " + (job_number + 1));
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
			/* ************************************************************************************************************************************************** */
			System.out.println("Reading line " + (line_count + 1) + ", the SEVENTH LINE of job: " + (job_number + 1));
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
			/* ************************************************************************************************************************************************** */
			line = lines.get(line_count); // Reads the divider line between jobs
			System.out.println("Reading line " + (line_count + 1));
			System.out.println(line + "\n");
			line_count++;
			/* ************************************************************************************************************************************************** */
		}

	private static void initialize_Job_Log()
		{
			try
				{
					log = new File(out_dir + "VIZ_LOG.txt");
					log_writer = new BufferedWriter(new FileWriter(log));
					log_writer.write("JED Cartesian Mode Visualization Driver version 1.0" + "\n");
					log_writer.write("Reference PDB file: " + PDB + "\n");
					log_writer.write("Eigenvector file: " + eigenvectors + "\n");
					log_writer.write("Modes file: " + modes + "\n");
					log_writer.write("Mode Maxes file: " + maxes + "\n");
					log_writer.write("Mode Mins file: " + mins + "\n");
					log_writer.write("Output directory: " + out_dir + "\n");
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
					log_writer.write("\nPerforming Cartesian Mode Visualization on Top  " + number_of_modes_viz + "  cPCA modes.\n");
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

	public static void main(String[] args)
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

					System.out.println("Performing Cartesian Mode Visualizations... ");
					startTime = System.nanoTime();

					get_Mode_Visualizations();
					write_VIZ_Log();

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
