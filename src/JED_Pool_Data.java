package jed;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.StringTokenizer;

import Jama.Matrix;

/**
 * JED class JED_Pool_Data: Constructs augmented matrices from the Coordinate Matrices from multiple trajectories.
 * 
 * The vector packing for n residues is by column: {X1...Xn,Y1...Yn,Z1...Zn}
 * 
 * The number of rows should be the same for all input matrices.
 * 
 * The number of columns in the output matrix is equal to the sum of all the input matrix columns.
 * 
 * Parameters are specified in the input file called "pool.txt", which is read when the program is run.
 * 
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
 * @author Dr. Charles David
 */

public class JED_Pool_Data
{
	static String line, input_path, out_dir, path;
	static int number_Of_Input_Lines, line_count, number_of_jobs, job_number, number_of_input_matrices, ROWS, COLS, Col_sum;
	static ArrayList<double[]> matrices;
	static ArrayList<String> lines, paths;
	static double[] in_matrix_vals;
	static int[] column_indices_start, column_indices_end;
	static Matrix augmented_matrix;
	static StringTokenizer sToken;
	static File Job_Log;
	static BufferedReader input_reader;
	static PrintWriter Job_Log_Writer;
	static boolean OK, check, isDir, exist;

	private static void read_input_file()
		{
			number_Of_Input_Lines = 0;
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
							System.err.println("THERE MUST BE AT LEAST 3 LINES OF PARAMETERS FOR ONE JOB.");
							System.err.println("Terminating program execution.");
							System.exit(0);
						}

				} catch (IOException e)
				{
					System.err.println("IOException thrown. Could not read the input file. Program will terminate.\n");
					e.printStackTrace();
					System.exit(0);
				}
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
			boolean OK = Test_Numeric_Type.test_Integer(test);
			if (OK) number_of_input_matrices = Integer.parseInt(test);
			System.out.println("\t# of Matrices to Combine = " + number_of_input_matrices);
			if (!OK)
				{
					System.err.println("The number of Matrices MUST be an integer.");
					System.err.println("Terminating program execution.");
					System.exit(0);
				}
			line_count++;
			/* ************************************************************************************************************************************************** */
			System.out.println("Reading line " + (line_count + 1) + ", the SECOND LINE of job: " + (job_number + 1));
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
			paths = new ArrayList<String>();
			for (int i = 0; i < number_of_input_matrices; i++)
				{
					System.out.println("Reading line " + (line_count + 1) + ", Matrix " + (i + 1) + " for job " + (job_number + 1));
					line = lines.get(line_count);
					sToken = new StringTokenizer(line);
					path = sToken.nextToken();
					paths.add(path);
					System.out.println("\tMatrix " + (i + 1) + " = " + path);
					exist = new File(path).exists();
					if (!exist)
						{
							System.err.println("The entered file does not exist: " + path);
							System.err.println("Terminating program execution.");
							System.exit(0);
						}
					line_count++;
				}
			/* ************************************************************************************************************************************************** */
			line = lines.get(line_count); // Reads the divider line between jobs
			System.out.println("Reading line " + (line_count + 1));
			System.out.println(line + "\n");
			line_count++;
			/* ************************************************************************************************************************************************** */
		}

	private static void initialize_Log()
		{
			Job_Log = new File(out_dir + "JED_Pool_Log_" + number_of_input_matrices + ".txt");
			try
				{
					Job_Log_Writer = new PrintWriter(new BufferedWriter(new FileWriter(Job_Log)));
				} catch (IOException e)
				{
					System.err.println("Could not create the Job Log file for Job Number " + job_number + 1 + ". Program terminating");
					e.printStackTrace();
					System.exit(0);
				}
		}

	private static void write_Log()
		{
			Job_Log_Writer.write("The number of matrices that were combined is: " + number_of_input_matrices + "\n");
			Job_Log_Writer.write("The output directory is: " + out_dir + "\n");
			Job_Log_Writer.write("The number of ROWS: " + ROWS + "\n");
			Job_Log_Writer.write("The total number of columns (frames): " + Col_sum + "\n");
			Job_Log_Writer.write("The pooled coordinate matrix was written to: " + path);
			Job_Log_Writer.flush();
			Job_Log_Writer.close();
		}

	/**
	 * Calls all private methods for the pooling of data.
	 * 
	 * @param args
	 *            Takes on argument: The full path to the input file.
	 */
	public static void main(String[] args)
		{
			System.out.println("Running the Pool Data Driver: ");
			System.out.println("Getting the input file: ");
			try
				{

					input_path = "pool.txt";
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
					check = new File(in_path).exists();
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

					matrices = new ArrayList<double[]>();
					column_indices_start = new int[number_of_input_matrices];
					column_indices_end = new int[number_of_input_matrices];
					Col_sum = 0;

					for (int i = 0; i < number_of_input_matrices; i++)
						{
							String input = paths.get(i);
							Matrix input_matrix = Matrix_IO.read_Matrix(input);
							ROWS = input_matrix.getRowDimension();
							COLS = input_matrix.getColumnDimension();
							column_indices_start[i] = Col_sum;
							column_indices_end[i] = Col_sum + COLS - 1;
							Col_sum += COLS;
							in_matrix_vals = input_matrix.getColumnPackedCopy();
							matrices.add(in_matrix_vals);
							in_matrix_vals = null;
							System.gc();
						}

					augmented_matrix = new Matrix(ROWS, Col_sum);
					for (int j = 0; j < number_of_input_matrices; j++)
						{
							int R1 = 0;
							int R2 = ROWS - 1;
							int C1 = column_indices_start[j];
							int C2 = column_indices_end[j];
							double[] mat = matrices.get(j);
							Matrix X = new Matrix(mat, ROWS);
							augmented_matrix.setMatrix(R1, R2, C1, C2, X);
							X = null;
							System.gc();
						}
					matrices.clear();

					path = out_dir + "Pooled_Coordinates_Matrix_" + number_of_input_matrices + ".txt";
					Matrix_IO.write_Matrix(augmented_matrix, path, 9, 3);

					augmented_matrix = null;
					System.gc();

					initialize_Log();
					write_Log();

					System.out.println("Job  " + (job_number + 1) + "  Done: ");
					System.out.println("-------------------------------------------------------------------------------------------------------------------------------------");
				}
		}
}
