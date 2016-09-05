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
 * JED class JED_Batch_Driver: Driver program for running JED in (sequential) batch mode.
 * Input file is "JED_Batch_Driver.txt"
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
 * along with this program. If not, see <http://www.gnu.org/license
 *
 * @author Dr. Charles David
 *
 */
public class JED_Batch_Driver
{

	static int number_Of_Input_Lines, line_count, num_of_jobs, job_num, number_of_modes_cartesian, number_of_modes_dist_pairs, number_of_residue_pairs, number_of_modes_viz,
			number_of_residues_cartesian, reference_column, ROWS, COLS;
	static String input_path, path, name, directory, out_dir, rl_cartesian, rl_dist_pairs, coord_input_file, description, pdb_ref_file, line, date, Q = "COV", R = "CORR", PC = "PCORR";
	static double trace_COV, trace_d_cov, trace_CORR, trace_d_corr, trace_PCORR, trace_d_pcorr, cond_cov, det_cov, rank_cov, cond_d_cov, det_d_cov, rank_d_cov, mode_amplitude, z_cut, percent;
	static List<Integer> residue_list, residue_list_dp1, residue_list_dp2, residue_list_dp_original1, residue_list_dp_original2, residues_read;
	static List<String> lines, pdb_file_names, chain_idents, chain_idents1, chain_idents2, chain_ids_read;
	static List<Double> transformed_conformation_rmsds, transformed_residue_rmsds, top_cartesian_eigenvalues_COV, top_distance_eigenvalues_COV, top_cartesian_eigenvalues_CORR,
			top_distance_eigenvalues_CORR, top_cartesian_eigenvalues_PCORR, top_distance_eigenvalues_PCORR;
	static double[] pca_mode_min_COV, pca_mode_max_COV, pca_mode_min_CORR, pca_mode_max_CORR, pca_mode_max_PCORR, pca_mode_min_PCORR, residue_distance_means, residue_distance_std_devs;
	static Matrix original_PDB_coordinates, subset_PDB_coordinates_cartesian, weighted_pca_modes_COV, square_pca_modes_COV, weighted_square_pca_modes_COV, transformed_PDB_coordinates_cartesian,
			distance_matrix, top_cartesian_evectors_COV, top_distance_evectors_COV, projections_COV, normed_projections_COV, weighted_pca_modes_CORR, square_pca_modes_CORR,
			weighted_square_pca_modes_CORR, top_cartesian_evectors_CORR, top_distance_evectors_CORR, projections_CORR, normed_projections_CORR, weighted_normed_projections_COV,
			weighted_normed_projections_CORR, pca_modes_CORR, pca_modes_COV, top_cartesian_evectors_PCORR, top_distance_evectors_PCORR, square_pca_modes_PCORR;
	static Vector<Atom> atoms;
	static File log;
	static BufferedReader input_reader;
	static BufferedWriter log_writer;
	static DateUtils now;
	static StringTokenizer sToken;
	static NumberFormat nf, nff;
	static RoundingMode rm;
	static boolean read_PDBs, do_cartesian, do_dist_pairs, do_mode_viz, do_multi, do_no_pca;
	static long startTime, endTime, totalTime;

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

				} catch (IOException e)
				{
					System.err.println("IOException thrown. Could not read the input file. Program will terminate.");
					e.printStackTrace();
					System.exit(0);
				}

		}

	private static void read_batch_parameters()
		{

			System.out.println("Reading line " + (line_count + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			String test = sToken.nextToken();
			boolean OK = Test_Numeric_Type.test_Integer(test);
			if (OK) num_of_jobs = Integer.parseInt(test);
			System.out.println("\tThe number of jobs =  " + num_of_jobs);
			if (!OK || num_of_jobs < 1)
				{
					System.err.println("Expected Number of Jobs to be a positive integer, but got: " + test);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			line_count++;
			/* ****************************************************************************************************************************************** */
			System.out.println("Reading line " + (line_count + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			int count = sToken.countTokens();
			if (count < 2)
				{
					System.err.println("Expected 2 parameters, $read and $multi, but only got: " + count);
					System.err.println(line);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}

			String read = sToken.nextToken();
			if (read.equals("1")) read_PDBs = true;
			System.out.println("\tREAD PDBs = " + read_PDBs);

			String multi_chain_analysis = sToken.nextToken();
			if (multi_chain_analysis.equals("1")) do_multi = true;
			System.out.println("\tMULTI CHAIN PDBs = " + do_multi);
			line_count++;
			if (read_PDBs)
				{
					percent = 0;
					z_cut = 0;
					do_no_pca = true;
					/* Reads the separator line between the batch parameters and the job parameters. */
					line = lines.get(line_count);
					System.out.println("Reading line " + (line_count + 1));
					System.out.println(line + "\n");
					line_count++;
					return;
				}
			/* ****************************************************************************************************************************************** */
			System.out.println("Reading line " + (line_count + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			count = sToken.countTokens();
			if (count < 2)
				{
					System.err.println("Expected 2 parameters, $percent and $Z_cutoff, but only got: " + count);
					System.err.println(line);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}

			test = sToken.nextToken();
			OK = Test_Numeric_Type.test_Double(test);
			if (OK) percent = Double.parseDouble(test);
			if (!OK || percent < 0 || percent >= 1)
				{
					System.err.println("Expected Percent to be a decimal in the range [0,1), but got: " + test);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			System.out.println("\tPercent = " + percent);
			test = sToken.nextToken();
			OK = Test_Numeric_Type.test_Double(test);
			if (OK) z_cut = Double.parseDouble(test);
			System.out.println("\tZ-Cutoff = " + z_cut);
			if (!OK || z_cut < 0)
				{
					System.err.println("Expected Z_cutoff to be zero or a positive decimal, but got: " + test);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			line_count++;
			/* ****************************************************************************************************************************************** */
			System.out.println("Reading line " + (line_count + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			count = sToken.countTokens();
			if (count < 4)
				{
					System.err.println("Expected AT LEAST 3 parameters: #Cartesian Modes, #Distance Pairs Modes, and #Modes Viz, but only got: " + count);
					System.err.println(line);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			test = sToken.nextToken();
			OK = Test_Numeric_Type.test_Integer(test);
			if (OK) number_of_modes_cartesian = Integer.parseInt(test);
			System.out.println("\t#Cartesian Modes = " + number_of_modes_cartesian);
			if (!OK)
				{
					System.err.println("The number of Cartesian Modes MUST be an integer.");
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			test = sToken.nextToken();
			OK = Test_Numeric_Type.test_Integer(test);
			if (OK) number_of_modes_dist_pairs = Integer.parseInt(test);
			System.out.println("\t#Residue Pair Modes = " + number_of_modes_dist_pairs);
			if (!OK)
				{
					System.err.println("The number of Residue Pair Modes MUST be an integer.");
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			test = sToken.nextToken();
			OK = Test_Numeric_Type.test_Integer(test);
			if (OK) number_of_modes_viz = Integer.parseInt(test);
			System.out.println("\t#Modes Viz = " + number_of_modes_viz);
			if (!OK)
				{
					System.err.println("The number of Modes to Visualize MUST be an integer.");
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			if (number_of_modes_viz > 0)
				{
					if (number_of_modes_viz > number_of_modes_cartesian)
						{
							System.err.println("The number of Cartesian Modes to Visualize is greater than the number of Cartesian Modes Available:");
							System.err.println("Number of Cartesian Modes to Visualize: " + number_of_modes_viz);
							System.err.println("Number of Cartesian Modes AVAILABLE: " + number_of_modes_cartesian);
							System.err.println("The Possible number of Cartesial Modes to Visualize is always <= Number of Cartesian Modes requested.");
							System.err.println("Setting the number of Cartesian Modes to Visualize to the number of Cartesian Modes Available:");
							number_of_modes_viz = number_of_modes_cartesian;
						}
					if (!sToken.hasMoreElements())
						{
							System.err.println("The number of Cartesian Modes to Visualize is > 0: " + number_of_modes_viz);
							System.err.println("Expected the Mode Amplitude parameter on this line (Field 4).");
							System.err.println("Setting the mode amplitude to the default value of 1.5");
							mode_amplitude = 1.5;
						}
					if (sToken.hasMoreElements())
						{
							test = sToken.nextToken();
							OK = Test_Numeric_Type.test_Double(test);
							if (OK) mode_amplitude = Double.parseDouble(test);
							System.out.println("\tMode Amplitude = " + mode_amplitude);
							if (!OK || mode_amplitude <= 0)
								{
									System.err.println("Mode Amplitude MUST be a decimal number > 0.");
									System.err.println("Suggested values range from 1.5 to 3.5 depending on the system");
									System.err.println("Setting mode apmlitude to default value of 1.5");
									mode_amplitude = 1.5;
								}
						}
				}
			line_count++;
			/* *********************************************************************************************************************************** */
			line = lines.get(line_count);
			System.out.println("Reading line " + (line_count + 1));
			System.out.println(line + "\n");
			line_count++;
			/* *********************************************************************************************************************************** */
			if (number_of_modes_cartesian > 0) do_cartesian = true;
			if (number_of_modes_dist_pairs > 0) do_dist_pairs = true;
			if (number_of_modes_viz > 0 & do_cartesian) do_mode_viz = true;
			/* *********************************************************************************************************************************** */
		}

	private static void initialize_Job_Log()
		{
			try
				{
					log = new File(out_dir + "JED_LOG.txt");
					log_writer = new BufferedWriter(new FileWriter(log));
					log_writer.write("Java Essential Dynamics" + "\n");
					log_writer.write("Results for job: " + (job_num + 1) + "\n");
					log_writer.write("Job Description: " + description + "\n");
					log_writer.write("Working directory: " + directory + "\n");
					log_writer.write("Output directory: " + out_dir + "\n");
					log_writer.write("READ PDBs =  " + read_PDBs + "\n");
					log_writer.write("MULTI CHAIN PDBs = " + do_multi + "\n");
					log_writer.write("Batch Parameters will apply to ALL jobs in this batch: " + "\n");
					log_writer.flush();

				} catch (IOException e)
				{
					System.err.println("IOException thrown. Could not write the JED_LOG file. Program will terminate.");
					e.printStackTrace();
					System.exit(0);
				}
		}

	private static void read_analytical_job_parameters()
		{

			if (line_count >= lines.size())
				{
					System.err.println("No more data in input file for remaining jobs in batch.");
					System.err.println("User specified too many jobs.");
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			System.out.println("Reading line " + (line_count + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			directory = sToken.nextToken();
			System.out.println("\tWorking Directory = " + directory);
			if (!(directory.endsWith("/") || directory.endsWith("\\")))
				{
					System.err.println("Expected the working Directory to end with / or \\\\, but got: " + line);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			boolean dir = new File(directory).isDirectory();
			if (!dir)
				{
					System.err.println("The entered directory does not exist as a directory: " + directory);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			line_count++;
			/* ****************************************************************************************************************************************** */
			System.out.println("Reading line " + (line_count + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			int count = sToken.countTokens();
			if (count < 2)
				{
					System.err.println("Expected 2 parameters: Job Description and Reference PDB File, but only got: " + count);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			description = sToken.nextToken();
			System.out.println("\tDescription = " + description);
			pdb_ref_file = sToken.nextToken();
			System.out.println("\tPDB Reference File = " + pdb_ref_file);
			boolean check = new File(directory + pdb_ref_file).exists();
			if (!check)
				{
					System.err.println("The entered PDB Reference File does not exist: " + directory + pdb_ref_file);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			out_dir = directory + "JED_RESULTS_" + description + "/";
			boolean exist = new File(out_dir).exists();
			if (!exist)
				{
					boolean success = (new File(out_dir)).mkdirs();
					if (!success)
						{
							System.out.println("Failed to create the output directory:  " + out_dir);
							System.exit(0);
						}
				}
			line_count++;
			/* ****************************************************************************************************************************************** */
			if (do_cartesian)
				{
					System.out.println("Reading line " + (line_count + 1));
					line = lines.get(line_count);
					sToken = new StringTokenizer(line);
					rl_cartesian = sToken.nextToken();
					System.out.println("\tCartesian Residue List = " + rl_cartesian);
					check = new File(directory + rl_cartesian).exists();
					if (!check)
						{
							System.err.println("The entered Cartesian Residue List File does not exist: " + directory + rl_cartesian);
							System.err.println("Terminating program execution.");
							System.exit(0);
						}
					line_count++;
				}
			/* ****************************************************************************************************************************************** */
			if (do_dist_pairs)
				{

					System.out.println("Reading line " + (line_count + 1));
					line = lines.get(line_count);

					sToken = new StringTokenizer(line);
					rl_dist_pairs = sToken.nextToken();
					System.out.println("\tDistance Residue Pair List = " + rl_dist_pairs);
					check = new File(directory + rl_dist_pairs).exists();
					if (!check)
						{
							System.err.println("The entered Distance Residue Pair List File does not exist: " + directory + rl_dist_pairs);
							System.err.println("Terminating program execution.");
							System.exit(0);
						}
					line_count++;
				}
			/* ****************************************************************************************************************************************** */
			if (!read_PDBs)
				{

					System.out.println("Reading line " + (line_count + 1));
					line = lines.get(line_count);

					sToken = new StringTokenizer(line);
					count = sToken.countTokens();
					if (count < 2)
						{
							System.err.println("Expected 2 parameters: Coordinate File Matrix and the Reference Column, but only got: " + count);
							System.err.println("Terminating program execution.");
							System.exit(0);
						}
					coord_input_file = sToken.nextToken();
					System.out.println("\tMatrix of Coordinates = " + coord_input_file);

					check = new File(directory + coord_input_file).exists();
					if (!check)
						{
							System.err.println("The entered Coordinates Matrix File does not exist: " + directory + coord_input_file);
							System.err.println("Terminating program execution.");
							System.exit(0);
						}
					String test = sToken.nextToken();
					boolean OK = Test_Numeric_Type.test_Integer(test);
					if (OK) reference_column = Integer.parseInt(test);
					System.out.println("\tReference Column = " + reference_column);
					if (!OK)
						{
							System.err.println("The Reference Column MUST be an integer.");
							System.err.println("Terminating program execution.");
							System.exit(0);
						}
					line_count++;
				}
			/* ****************************************************************************************************************************************** */
			line = lines.get(line_count);
			System.out.println("Reading line " + (line_count + 1));
			System.out.println(line + "\n");
			line_count++;
			/* ****************************************************************************************************************************************** */
			initialize_Job_Log();
		}

	private static void read_preprocessing_job_parameters()
		{

			if (line_count >= lines.size())
				{
					System.err.println("No more data in input file for remaining jobs in batch.");
					System.err.println("User specified too many jobs.");
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			System.out.println("Reading line " + (line_count + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			directory = sToken.nextToken();
			System.out.println("\tWorking Directory = " + directory);
			if (!(directory.endsWith("/") || directory.endsWith("\\")))
				{
					System.err.println("Expected the working Directory to end with / or \\\\, but got: " + line);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			boolean dir = new File(directory).isDirectory();
			if (!dir)
				{
					System.err.println("The entered directory does not exist as a directory: " + directory);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			line_count++;
			/* ****************************************************************************************************************************************** */
			System.out.println("Reading line " + (line_count + 1));
			line = lines.get(line_count);
			sToken = new StringTokenizer(line);
			int count = sToken.countTokens();
			if (count < 2)
				{
					System.err.println("Expected 2 parameters: Job Description and Reference PDB File, but only got: " + count);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			description = sToken.nextToken();
			System.out.println("\tDescription = " + description);
			pdb_ref_file = sToken.nextToken();
			System.out.println("\tPDB Reference File = " + pdb_ref_file);
			boolean check = new File(directory + pdb_ref_file).exists();
			if (!check)
				{
					System.err.println("The entered PDB Reference File does not exist: " + directory + pdb_ref_file);
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			out_dir = directory + "JED_RESULTS_" + description + "/";
			boolean exist = new File(out_dir).exists();
			if (!exist)
				{
					boolean success = (new File(out_dir)).mkdirs();
					if (!success)
						{
							System.out.println("Failed to create the output directory:  " + out_dir);
							System.exit(0);
						}
				}
			line_count++;
			/* ****************************************************************************************************************************************** */
			line = lines.get(line_count);
			System.out.println("Reading line " + (line_count + 1));
			System.out.println(line + "\n");
			line_count++;
			/* ****************************************************************************************************************************************** */
			initialize_Job_Log();
		}
	/* ****************************************************************************************************************************************** */

	private static void read_pdb_files()
		{
			JED_Read_PDBs rPDBs = new JED_Read_PDBs(directory, description, pdb_ref_file);

			rPDBs.read_PDBs();

			original_PDB_coordinates = rPDBs.get_Original_PDB_Coordinates();
			pdb_file_names = rPDBs.get_PDB_File_Names();
			reference_column = rPDBs.get_Reference_Column();
			residues_read = rPDBs.get_Residues_Read();
			ROWS = original_PDB_coordinates.getRowDimension();
			COLS = original_PDB_coordinates.getColumnDimension();
			if (!do_multi)
				{
					name = "All_PDB_Residues_JED.txt";
					path = out_dir + name;
					List_IO.write_Integer_List(residues_read, path);
					write_Read_Log();
				}
			if (do_multi)
				{
					chain_ids_read = rPDBs.get_Chain_IDs_Read();
					name = "All_PDB_Residues_Multi_JED.txt";
					path = out_dir + name;
					List_IO.write_String_Integer_List(chain_ids_read, residues_read, path);
					write_Read_Log();
				}
			System.gc();
		}

	private static void write_Read_Log()
		{
			try
				{
					log_writer.write("\nPerforming Preliminary Processing Run.\n");
					log_writer.write("\nReading alpha carbon coordinates from all PDB files in the Working Directory" + "\n");
					log_writer.write("The number of PDB files read = " + COLS + "\n");
					log_writer.write("The number of residues found in the PDB files = " + (ROWS / 3) + "\n");
					log_writer.write("Coordinates Matrix created: 'original_PDB_Coordinates'" + "\n");
					log_writer.write("The dimension of the coordinates matrix is: " + ROWS + "  by  " + COLS + "\n");
					log_writer.write("The transformed PDB coordinates were obtained by quaternion least-squares alignment to the reference structure.\n");
					log_writer.write("PDB reference structure is " + pdb_ref_file + "\n");
					log_writer.write("The column in the Coordinates Matrix that corresponds to the entered reference structure is:\n");
					log_writer.write("\tREFERENCE COLUMN =  " + reference_column + "\n\n");
					if (!do_multi) log_writer.write("The file '" + name + "' contains all residue numbers found in the PDB file, formatted for Single Chain Analysis Input." + "\n");
					if (do_multi) log_writer.write("The file '" + name + "' contains all chainID-residue number pairs found in the PDB file, formatted for Multi Chain Analysis Input." + "\n");
					log_writer.flush();

				} catch (IOException e)
				{
					System.err.println("Could not write the JED_LOG file: " + out_dir + "JED_LOG.txt");
					System.err.println("Program terminating.\n");
					e.printStackTrace();
					System.exit(0);
				}
			System.gc();
		}

	private static void read_coordinate_file()
		{

			JED_Get_Coordinates_from_Matrix mat_coords = new JED_Get_Coordinates_from_Matrix(directory, coord_input_file);
				{
					original_PDB_coordinates = mat_coords.get_Original_PDB_coordinates();
					ROWS = original_PDB_coordinates.getRowDimension();
					COLS = original_PDB_coordinates.getColumnDimension();
				}

				{
					try
						{
							log_writer.write("\nThe alpha carbon coordinates were obtained from coordinates matrix file: " + coord_input_file + "\n");
							log_writer.write("The dimension of the coordinates matrix is = " + ROWS + " by " + COLS + "\n");
							log_writer.write("Total number of residues in matrix = " + (ROWS / 3) + "\n");
							log_writer.write("Total number of conformations in matrix = " + COLS + "\n");
							log_writer.write("Transformed PDB coordinates obtained by quaternion least-squares alignment to the reference structure." + "\n");
							log_writer.write("PDB reference structure is: " + pdb_ref_file + "\n");
							log_writer.write("Reference Column in matrix = " + reference_column + "\n");
							log_writer.flush();

						} catch (IOException e)
						{
							System.err.println("Could not write the JED_LOG file: " + out_dir + "JED_LOG.txt");
							System.err.println("Program terminating.\n");
							e.printStackTrace();
							System.exit(0);
						}
					System.gc();
				}
		}

	private static void do_cPCA()
		{

			JED_Do_Cartesian cSS = new JED_Do_Cartesian(directory, description, pdb_ref_file, rl_cartesian, reference_column, number_of_modes_cartesian, original_PDB_coordinates);

			cSS.set_percent_cutoff(percent);
			cSS.set_z_cutoff(z_cut);

			if (!do_multi) cSS.do_Cartesian();
			if (do_multi) cSS.do_Cartesian_Multi();

			number_of_residues_cartesian = cSS.getNumber_of_residues();
			atoms = cSS.get_atoms();
			trace_COV = cSS.get_Trace_COV();
			trace_CORR = cSS.get_Trace_CORR();
			trace_PCORR = cSS.get_Trace_PCORR();
			cond_cov = cSS.get_cond_COV();
			det_cov = cSS.get_Det_COV();
			rank_cov = cSS.get_Rank_COV();
			top_cartesian_eigenvalues_COV = cSS.getTop_cartesian_eigenvalues_COV();
			top_cartesian_eigenvalues_CORR = cSS.getTop_cartesian_eigenvalues_CORR();
			top_cartesian_eigenvalues_PCORR = cSS.getTop_cartesian_eigenvalues_PCORR();
			top_cartesian_evectors_COV = cSS.getTop_cartesian_evectors_COV();
			top_cartesian_evectors_CORR = cSS.getTop_cartesian_evectors_CORR();
			top_cartesian_evectors_PCORR = cSS.getTop_cartesian_evectors_PCORR();
			square_pca_modes_COV = cSS.getSquare_pca_modes_COV();
			square_pca_modes_CORR = cSS.getSquare_pca_modes_CORR();
			square_pca_modes_PCORR = cSS.getSquare_pca_modes_PCORR();
			pca_mode_max_COV = cSS.getPca_mode_max_COV();
			pca_mode_max_CORR = cSS.getPca_mode_max_CORR();
			pca_mode_max_PCORR = cSS.getPca_mode_max_PCORR();
			pca_mode_min_COV = cSS.getPca_mode_min_COV();
			pca_mode_min_CORR = cSS.getPca_mode_min_CORR();
			pca_mode_min_PCORR = cSS.getPca_mode_min_PCORR();

			System.gc();

			write_cPCA_Log();
			get_cPCA_SSA();
		}

	private static void get_cPCA_SSA()
		{
			JED_Get_Subspace_Analysis jssa1 = new JED_Get_Subspace_Analysis(directory, description, "cPCA", top_cartesian_evectors_COV, top_cartesian_evectors_CORR, "COV_vs_CORR");
				{
					jssa1.get_SSA_JED();
					jssa1.get_FSSA_Iterated_JED();
					System.gc();
				}
			JED_Get_Subspace_Analysis jssa2 = new JED_Get_Subspace_Analysis(directory, description, "cPCA", top_cartesian_evectors_COV, top_cartesian_evectors_PCORR, "COV_vs_PCORR");
				{
					jssa2.get_SSA_JED();
					jssa2.get_FSSA_Iterated_JED();
					System.gc();
				}
			JED_Get_Subspace_Analysis jssa3 = new JED_Get_Subspace_Analysis(directory, description, "cPCA", top_cartesian_evectors_CORR, top_cartesian_evectors_PCORR, "CORR_vs_PCORR");
				{
					jssa3.get_SSA_JED();
					jssa3.get_FSSA_Iterated_JED();
					System.gc();
				}
		}

	private static void write_cPCA_Log()
		{
			try
				{
					log_writer.write("\nPERFORMING cPCA, Computing Top " + number_of_modes_cartesian + " modes." + "\n\n");
					log_writer.write("Residue list for Cartesian subset:  " + rl_cartesian + "\n");
					log_writer.write("Number of residues in Cartesian subset: " + number_of_residues_cartesian + "\n");
					if (percent > 0) log_writer.write("The transformed data was trimmed by removing " + (percent * 100) + " percent of the samples based on conformation RMSD." + "\n");
					else log_writer.write("No samples were removed from the data" + "\n");
					if (z_cut > 0) log_writer.write("The coordinate values with Z-scores beyond |" + z_cut + "| were set to their mean value.\n\n");
					else log_writer.write("No coordinate outliers were adjusted.\n");
					if (trace_COV >= 1.00) log_writer.write("Trace of the Covariance Matrix = " + nff.format(trace_COV) + "\n");
					if (trace_COV < 1.00) log_writer.write("Trace of the Covariance Matrix = " + nf.format(trace_COV) + "\n");
					log_writer.write("Condition Number of the Covariance Matrix = " + nff.format(cond_cov) + "\n");
					log_writer.write("Determinant of the Covariance Matrix = " + (det_cov) + "\n");
					log_writer.write("Rank of the Covariance Matrix = " + nff.format(rank_cov) + "\n");
					log_writer.write("Trace of the Correlation Matrix = " + nff.format(trace_CORR) + "\n");
					log_writer.write("Trace of the Partial Correlation Matrix = " + nff.format(trace_PCORR) + "\n");
					log_writer.write("PDB file with B-factors replaced by residue RMSDs: ss_" + number_of_residues_cartesian + "_RMSF_edited.pdb" + "\n");
					log_writer.write("The DVPs (PCs) from the 3 different models were calculated using:" + "\n");
					log_writer.write("Standard dot product(dp), normed dp, weighted dp (by eigenvalue), and weighted normed dp" + "\n");
					if (number_of_modes_cartesian < 2) log_writer.write("Since the number of modes is less than 2, no FE calculation will be done." + "\n");
					if (number_of_modes_cartesian >= 2)
						{
							log_writer.write("The Free energy (delta_G) was determined using the first 2 standard DVPs as order parameters." + "\n");
							log_writer.write("Probabilities were calculated using 2D Kernel Density Estimation with Gaussian kernels." + "\n");
						}
					log_writer.write("Subspace analysis was done comparing the top vector spaces from the 3 different models." + "\n");
					log_writer.write("Comparators include RMSIP and Principle Angles, for the essential subspace and iterated comparisons from dim 1 to " + number_of_modes_cartesian + "\n");
					log_writer.write("Additional log files can be found in the /SSA directory tree." + "\n\n");
					log_writer.flush();

				} catch (IOException e)
				{
					System.err.println("Could not write the JED_LOG file: " + out_dir + "JED_LOG.txt");
					System.err.println("Program terminating.\n");
					e.printStackTrace();
					System.exit(0);
				}
			System.gc();
		}

	private static void do_dpPCA()
		{
			JED_Do_Dist_Pairs dpSS = new JED_Do_Dist_Pairs(directory, description, pdb_ref_file, rl_dist_pairs, reference_column, number_of_modes_dist_pairs, original_PDB_coordinates);
			dpSS.set_z_cutoff(z_cut);

			if (!do_multi) dpSS.do_Dist();
			if (do_multi) dpSS.do_Dist_Multi();

			number_of_residue_pairs = dpSS.getNumber_of_pairs();
			residue_distance_means = dpSS.getResidue_distance_means();
			residue_distance_std_devs = dpSS.getResidue_distance_std_devs();
			residue_list_dp1 = dpSS.getResidue_list_dist1();
			residue_list_dp2 = dpSS.getResidue_list_dist2();
			residue_list_dp_original1 = dpSS.getResidue_list_dist_orig1();
			residue_list_dp_original2 = dpSS.getResidue_list_dist_orig2();
			chain_idents1 = dpSS.getChain_idents1();
			chain_idents2 = dpSS.getChain_idents2();

			trace_d_cov = dpSS.getTrace_dist_COV();
			cond_d_cov = dpSS.getCond_cov();
			det_d_cov = dpSS.getDet_cov();
			rank_d_cov = dpSS.getRank_cov();
			trace_d_corr = dpSS.getTrace_dist_CORR();
			trace_d_pcorr = dpSS.getTrace_dist_PCORR();

			top_distance_evectors_COV = dpSS.getTop_distance_evectors_COV();
			top_distance_evectors_CORR = dpSS.getTop_distance_evectors_CORR();
			top_distance_evectors_PCORR = dpSS.getTop_distance_evectors_PCORR();

			System.gc();

			if (!do_multi) write_dpPCA_Log_Single();
			if (do_multi) write_dpPCA_Log_Multi();

			get_dpPCA_SSA();
		}

	private static void get_dpPCA_SSA()
		{
			JED_Get_Subspace_Analysis jssa1 = new JED_Get_Subspace_Analysis(directory, description, "dpPCA", top_cartesian_evectors_COV, top_cartesian_evectors_CORR, "COV_vs_CORR");
				{
					jssa1.get_SSA_JED();
					jssa1.get_FSSA_Iterated_JED();
					System.gc();
				}
			JED_Get_Subspace_Analysis jssa2 = new JED_Get_Subspace_Analysis(directory, description, "dpPCA", top_cartesian_evectors_COV, top_cartesian_evectors_PCORR, "COV_vs_PCORR");
				{
					jssa2.get_SSA_JED();
					jssa2.get_FSSA_Iterated_JED();
					System.gc();
				}
			JED_Get_Subspace_Analysis jssa3 = new JED_Get_Subspace_Analysis(directory, description, "dpPCA", top_cartesian_evectors_CORR, top_cartesian_evectors_PCORR, "CORR_vs_PCORR");
				{
					jssa3.get_SSA_JED();
					jssa3.get_FSSA_Iterated_JED();
					System.gc();
				}
		}

	private static void write_dpPCA_Log_Multi()
		{
			try
				{
					log_writer.write("\nPERFORMING dpPCA, Computing Top " + number_of_modes_dist_pairs + " modes.\n");
					log_writer.write("\nResidue Pair list:  " + rl_dist_pairs + "\n");
					log_writer.write("Number of residues pairs: " + number_of_residue_pairs + "\n");
					if (z_cut > 0) log_writer.write("The distance values with Z-scores beyond |" + z_cut + "| were set to their mean value.\n\n");
					else log_writer.write("No coordinate outliers were adjusted.\n");
					if (trace_d_cov >= 1.00) log_writer.write("Trace of the Covariance Matrix = " + nff.format(trace_d_cov) + "\n");
					if (trace_d_cov < 1.00) log_writer.write("Trace of the Covariance Matrix = " + nf.format(trace_d_cov) + "\n");
					log_writer.write("Condition Number of the Covariance Matrix = " + nff.format(cond_d_cov) + "\n");
					log_writer.write("Determinant of the Covariance Matrix = " + nff.format(det_d_cov) + "\n");
					log_writer.write("Rank of the Covariance Matrix = " + nff.format(rank_d_cov) + "\n");
					log_writer.write("Trace of the Correlation Matrix = " + nff.format(trace_d_corr) + "\n");
					log_writer.write("Trace of the Partial Correlation Matrix = " + nff.format(trace_d_pcorr) + "\n");
					log_writer.write("\nMEANs and STANDARD DEVIATIONs for the Residue Pair Distances: " + "\n\n");
					log_writer.write(String.format("%-12s%-16s%-16s%-16s%n", "Res1", "Res2", "Mean", "Std_Dev"));
					for (int i = 0; i < number_of_residue_pairs; i++)
						{
							log_writer.write(String.format("%-12s%-16s%-16s%-16s%n", chain_idents1.get(i) + residue_list_dp_original1.get(i), chain_idents2.get(i) + residue_list_dp_original2.get(i),
									nf.format(residue_distance_means[i]), nf.format(residue_distance_std_devs[i])));
						}
					log_writer.write("\nThe DVPs (PCs) from from the 3 different models were calculated using:" + "\n");
					log_writer.write("Standard dot product(dp), normed dp, weighted dp (by eigenvalue), and weighted normed dp." + "\n");
					if (number_of_modes_dist_pairs < 2) log_writer.write("Since the number of modes is less than 2, no FE calculation will be done." + "\n");
					if (number_of_modes_dist_pairs >= 2)
						{
							log_writer.write("The Free energy (delta_G) was determined using the first 2 standard DVPs as order parameters." + "\n");
							log_writer.write("Probabilities were calculated using 2D Kernel Density Estimation with Gaussian kernels." + "\n");
						}
					log_writer.write("Subspace analysis was done to compare the top vector spaces from the 3 different models." + "\n");
					log_writer.write("Comparators include RMSIP and Principle Angles, for the essential subspace and iterated comparisons from dim 1 to " + number_of_modes_dist_pairs + "\n");
					log_writer.write("Additional log files can be found in the /SSA directory tree." + "\n\n");
					log_writer.flush();

				} catch (IOException e)
				{
					System.err.println("Could not write the JED_LOG file: " + out_dir + "JED_LOG.txt");
					System.err.println("Program terminating.\n");
					e.printStackTrace();
					System.exit(0);
				}
			System.gc();
		}

	private static void write_dpPCA_Log_Single()
		{
			try
				{
					log_writer.write("\nPERFORMING dpPCA, Computing Top " + number_of_modes_dist_pairs + " modes.\n");
					if (number_of_modes_dist_pairs < 2) log_writer.write("\nSince the number of modes is less than 2, no FE calculation will be done." + "\n");
					log_writer.write("\nResidue Pair list:  " + rl_dist_pairs + "\n");
					log_writer.write("Number of residue pairs: " + number_of_residue_pairs + "\n");
					if (z_cut > 0) log_writer.write("The distance values with Z-scores beyond |" + z_cut + "| were set to their mean value.\n\n");
					else log_writer.write("No coordinate outliers were adjusted.\n");
					if (trace_d_cov >= 1.00) log_writer.write("Trace of the Covariance Matrix = " + nff.format(trace_d_cov) + "\n");
					if (trace_d_cov < 1.00) log_writer.write("Trace of the Covariance Matrix = " + nf.format(trace_d_cov) + "\n");
					log_writer.write("Condition Number of the Covariance Matrix = " + nff.format(cond_d_cov) + "\n");
					log_writer.write("Determinant of the Covariance Matrix = " + nff.format(det_d_cov) + "\n");
					log_writer.write("Rank of the Covariance Matrix = " + nff.format(rank_d_cov) + "\n");
					log_writer.write("Trace of the Correlation Matrix = " + nff.format(trace_d_corr) + "\n");
					log_writer.write("Trace of the Partial Correlation Matrix = " + nff.format(trace_d_pcorr) + "\n");
					log_writer.write("\nMEANs and STANDARD DEVIATIONs for the Residue Pair Distances: " + "\n\n");
					log_writer.write(String.format("%-12s%-16s%-16s%-16s%n", "Res1", "Res2", "Mean", "Std_Dev"));
					for (int i = 0; i < number_of_residue_pairs; i++)
						{
							log_writer.write(String.format("%-12s%-16s%-16s%-16s%n", residue_list_dp_original1.get(i), residue_list_dp_original2.get(i), nf.format(residue_distance_means[i]),
									nf.format(residue_distance_std_devs[i])));
						}
					log_writer.write("\nThe DVPs (PCs) from from the 3 different models were calculated using:" + "\n");
					log_writer.write("Standard dot product(dp), normed dp, weighted dp (by eigenvalue), and weighted normed dp." + "\n");
					if (number_of_modes_dist_pairs >= 2)
						{
							log_writer.write("The Free energy (delta_G) was determined using the first 2 standard DVPs as order parameters." + "\n");
							log_writer.write("Probabilities were calculated using 2D Kernel Density Estimation with Gaussian kernels." + "\n");
						}
					log_writer.write("Subspace analysis was done to compare the top vector spaces from the 3 different models." + "\n");
					log_writer.write("Comparators include RMSIP and Principle Angles, for the essential subspace and iterated comparisons from dim 1 to " + number_of_modes_dist_pairs + "\n");
					log_writer.write("Additional log files can be found in the /SSA directory tree." + "\n\n");
					log_writer.flush();

				} catch (IOException e)
				{
					System.err.println("Could not write the JED_LOG file: " + out_dir + "JED_LOG.txt");
					System.err.println("Program terminating.\n");
					e.printStackTrace();
					System.exit(0);
				}
			System.gc();
		}

	private static void do_NO_PCA()
		{
			if (!do_multi)
				{
					JED_Do_No_PCA npca = new JED_Do_No_PCA(directory, rl_cartesian, residues_read, description, pdb_ref_file, original_PDB_coordinates, reference_column);
						{
							npca.set_percent_cutoff(percent);
							npca.set_z_cutoff(z_cut);
							npca.do_No_PCA();
							transformed_PDB_coordinates_cartesian = npca.getTransformed_subset_PDB_coordinates();
							transformed_conformation_rmsds = npca.get_conformation_rmsds();
							transformed_residue_rmsds = npca.get_residue_rmsd_list();
							atoms = npca.getAtoms();
							number_of_residues_cartesian = npca.getNumber_of_residues();

							original_PDB_coordinates = null;
							transformed_PDB_coordinates_cartesian = null;

							System.gc();

							write_NO_PCA_Log();
						}
				}
			if (do_multi)
				{
					JED_Do_No_PCA npca = new JED_Do_No_PCA(directory, rl_cartesian, chain_ids_read, residues_read, description, pdb_ref_file, original_PDB_coordinates, reference_column);
						{
							npca.set_percent_cutoff(percent);
							npca.set_z_cutoff(z_cut);
							npca.do_No_PCA_Multi();
							transformed_PDB_coordinates_cartesian = npca.getTransformed_subset_PDB_coordinates();
							transformed_conformation_rmsds = npca.get_conformation_rmsds();
							transformed_residue_rmsds = npca.get_residue_rmsd_list();
							atoms = npca.getAtoms();
							number_of_residues_cartesian = npca.getNumber_of_residues();

							original_PDB_coordinates = null;
							transformed_PDB_coordinates_cartesian = null;

							System.gc();

							write_NO_PCA_Log();
						}
				}
		}

	private static void write_NO_PCA_Log()
		{
			try
				{
					log_writer.write("\nConformation RMSDs and Residue RMSDs (RMSF) were calculated." + "\n");
					log_writer.write("\nPDB file with B-factors replaced by residue RMSDs: ss_" + number_of_residues_cartesian + "_RMSF_edited.pdb" + "\n");
					log_writer.write("\nVariable Z-scores were calculated.\n");

					log_writer.flush();

				} catch (IOException e)
				{
					System.err.println("Could not write the JED_LOG file: " + out_dir + "JED_LOG.txt");
					System.err.println("Program terminating.\n");
					e.printStackTrace();
					System.exit(0);
				}
			System.gc();
		}

	private static void do_Mode_Visualization()
		{
			JED_PCA_Mode_Vizualization cov = new JED_PCA_Mode_Vizualization(directory, description, atoms, top_cartesian_eigenvalues_COV, top_cartesian_evectors_COV, square_pca_modes_COV,
					pca_mode_max_COV, pca_mode_min_COV, mode_amplitude, number_of_modes_viz, Q);
			cov.get_Mode_Visualizations_SS();
			cov.get_Essential_Visualization();

			System.gc();

			JED_PCA_Mode_Vizualization corr = new JED_PCA_Mode_Vizualization(directory, description, atoms, top_cartesian_eigenvalues_CORR, top_cartesian_evectors_CORR, square_pca_modes_CORR,
					pca_mode_max_CORR, pca_mode_min_CORR, mode_amplitude, number_of_modes_viz, R);
			corr.get_Mode_Visualizations_SS();
			corr.get_Essential_Visualization();

			JED_PCA_Mode_Vizualization pcorr = new JED_PCA_Mode_Vizualization(directory, description, atoms, top_cartesian_eigenvalues_PCORR, top_cartesian_evectors_PCORR, square_pca_modes_PCORR,
					pca_mode_max_PCORR, pca_mode_min_PCORR, mode_amplitude, number_of_modes_viz, PC);
			pcorr.get_Mode_Visualizations_SS();
			pcorr.get_Essential_Visualization();

			System.gc();

			write_VIZ_Log();
		}

	private static void write_VIZ_Log()
		{
			try
				{
					log_writer.write("Performing Cartesian Mode Visualization on Top  " + number_of_modes_viz + "  cPCA modes.\n");
					log_writer.write("Sets of 20 structures were generated to animate each selected cPCA mode, for the COV, CORR, and PCORR PCA models.\n");
					log_writer.write("Atoms of each residue were perturbed along the mode eigenvector using a sine function ranging from 0 to 2PI.\n");
					log_writer.write("A Pymol(TM) script was generated for each mode to play the mode structures as a movie.\n");
					log_writer.write("MODE AMPLITUDE = " + nf.format(mode_amplitude) + "\n");
					log_writer.flush();

				} catch (IOException e)
				{
					System.err.println("Could not write the JED_LOG file: " + out_dir + "JED_LOG.txt");
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
			nff = NumberFormat.getInstance();
			nff.setRoundingMode(rm);
			nff.setMaximumFractionDigits(0);
			nff.setMinimumFractionDigits(0);

			System.out.println("JED: Java Essential Dynamics version 1.0 Batch Driver: ");
			System.out.println("Getting the input file: ");
			try
				{

					input_path = "JED_Batch_Driver.txt";
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
							System.err.println("The entered Input File could not be found: " + in_path);
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

			System.out.println("Reading Batch Parameters: ");
			read_batch_parameters();

			for (job_num = 0; job_num < num_of_jobs; job_num++)
				{
					System.out.println("Now Running Job Number " + (job_num + 1));

					if (read_PDBs)
						{
							read_preprocessing_job_parameters();
							System.out.println("Reading PDB files... ");
							startTime = System.nanoTime();
							read_pdb_files();
							endTime = System.nanoTime();
							totalTime = endTime - startTime;
							System.out.println("Done. (" + nf.format(totalTime / 1000000000.0) + " seconds)");
						}

					if (!read_PDBs)
						{
							read_analytical_job_parameters();
							System.out.println("Reading Coordinate File... ");
							startTime = System.nanoTime();
							read_coordinate_file();
							endTime = System.nanoTime();
							totalTime = endTime - startTime;
							System.out.println("Done. (" + nf.format(totalTime / 1000000000.0) + " seconds)");
						}

					if (do_cartesian)
						{
							System.out.println("Performing cPCA analysis... ");
							startTime = System.nanoTime();
							do_cPCA();
							endTime = System.nanoTime();
							totalTime = endTime - startTime;
							System.out.println("Done. (" + nf.format(totalTime / 1000000000.0) + " seconds)");
						}
					if (do_dist_pairs)
						{
							System.out.println("Performing dpPCA analysis... ");
							startTime = System.nanoTime();
							do_dpPCA();
							endTime = System.nanoTime();
							totalTime = endTime - startTime;
							System.out.println("Done. (" + nf.format(totalTime / 1000000000.0) + " seconds)");
						}
					if (do_no_pca)
						{
							System.out.println("Pre-Processing PDB files for JED Analysis... ");
							startTime = System.nanoTime();
							do_NO_PCA();
							endTime = System.nanoTime();
							totalTime = endTime - startTime;
							System.out.println("Done. (" + nf.format(totalTime / 1000000000.0) + " seconds)");
						}
					if (do_mode_viz)
						{
							System.out.println("Performing mode visualization... ");
							startTime = System.nanoTime();
							do_Mode_Visualization();
							endTime = System.nanoTime();
							totalTime = endTime - startTime;
							System.out.println("Done. (" + nf.format(totalTime / 1000000000.0) + " seconds)");
						}
					date = DateUtils.now();
					try
						{
							log_writer.write("\nAnalysis completed: " + date);
							log_writer.close();

						} catch (IOException e)
						{
							System.err.println("IOException thrown. Could not write the JED_LOG file. Program will terminate.");
							e.printStackTrace();
							System.exit(0);
						}

					System.out.println("JED Job " + (job_num + 1) + " completed: " + date);
					System.out.println("--------------------------------------------------------------------------------------------------------------------------------");
				}
			date = DateUtils.now();
			System.out.println("JED batch of " + num_of_jobs + " jobs completed: " + date + "\n");
		}
}
