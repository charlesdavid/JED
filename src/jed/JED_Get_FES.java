package jed;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.Arrays;

import Jama.Matrix;

public class JED_Get_FES
{
	String Q = "COV", R = "CORR", PC = "PCORR", model, type;
	long startTime, endTime, totalTime;
	static String directory, out_directory, description, file_name_head, path, date;
	static int ROWS, COLS;
	static Integer OP1, OP2, number_of_points, OP_offset, KDE_offset;
	static double cellsize;
	static double[] KDE_bounds, KDE_bandwidths;
	static Matrix FE, delta_vector_series, order_param1, order_param2;
	static File LOG;
	static PrintWriter LOG_writer;
	static NumberFormat nf;
	boolean check, isDir, exist, success;
	KernelDensityEstimate2d KDE;

	public JED_Get_FES(Matrix DVPs, int op1, int op2, int num_points, int offset, double size, String dir, String desc, String PCA_type, String PCA_model)
		{
			delta_vector_series = DVPs;
			OP1 = op1;
			OP2 = op2;
			number_of_points = num_points;
			OP_offset = offset;
			cellsize = size;
			directory = dir;
			description = desc;
			type = PCA_type;
			model = PCA_model;

			ROWS = delta_vector_series.getRowDimension();
			COLS = delta_vector_series.getColumnDimension();
			int row_index1 = OP_offset;
			int row_index2 = OP_offset + number_of_points - 1;
			if (row_index2 > ROWS - 1) row_index2 = ROWS - 1;
			order_param1 = delta_vector_series.getMatrix(row_index1, row_index2, OP1, OP1);
			order_param2 = delta_vector_series.getMatrix(row_index1, row_index2, OP2, OP2);

			nf = NumberFormat.getInstance();
			nf.setMaximumFractionDigits(3);
			nf.setMinimumFractionDigits(3);
			nf.setRoundingMode(RoundingMode.HALF_UP);

			if (type.equals("cPCA"))
				{
					if (model.equals(Q))
						{
							out_directory = directory + "JED_RESULTS_" + description + File.separatorChar + "cPCA" + File.separatorChar + "COV" + File.separatorChar + "FES" + File.separatorChar;
							exist = new File(out_directory).exists();
							if (!exist)
								{
									success = new File(out_directory).mkdirs();
									if (!success) System.err.println("Could not create the output directory: " + out_directory);
								}
						}
					if (model.equals(R))
						{
							out_directory = directory + "JED_RESULTS_" + description + File.separatorChar + "cPCA" + File.separatorChar + "CORR" + File.separatorChar + "FES" + File.separatorChar;
							exist = new File(out_directory).exists();
							if (!exist)
								{
									success = new File(out_directory).mkdirs();
									if (!success) System.err.println("Could not create the output directory: " + out_directory);
								}
						}
					if (model.equals(PC))
						{
							out_directory = directory + "JED_RESULTS_" + description + File.separatorChar + "cPCA" + File.separatorChar + "PCORR" + File.separatorChar + "FES" + File.separatorChar;
							exist = new File(out_directory).exists();
							if (!exist)
								{
									success = new File(out_directory).mkdirs();
									if (!success) System.err.println("Could not create the output directory: " + out_directory);
								}
						}
				}
			if (type.equals("dpPCA"))
				{
					if (model.equals(Q))
						{
							out_directory = directory + "JED_RESULTS_" + description + File.separatorChar + "dpPCA" + File.separatorChar + "COV" + File.separatorChar + "FES" + File.separatorChar;
							exist = new File(out_directory).exists();
							if (!exist)
								{
									success = new File(out_directory).mkdirs();
									if (!success) System.err.println("Could not create the output directory: " + out_directory);
								}
						}
					if (model.equals(R))
						{
							out_directory = directory + "JED_RESULTS_" + description + File.separatorChar + "dpPCA" + File.separatorChar + "CORR" + File.separatorChar + "FES" + File.separatorChar;
							exist = new File(out_directory).exists();
							if (!exist)
								{
									success = new File(out_directory).mkdirs();
									if (!success) System.err.println("Could not create the output directory: " + out_directory);
								}
						}
					if (model.equals(PC))
						{
							out_directory = directory + "JED_RESULTS_" + description + File.separatorChar + "dpPCA" + File.separatorChar + "PCORR" + File.separatorChar + "FES" + File.separatorChar;
							exist = new File(out_directory).exists();
							if (!exist)
								{
									success = new File(out_directory).mkdirs();
									if (!success) System.err.println("Could not create the output directory: " + out_directory);
								}
						}
				}
			file_name_head = out_directory + description + "_FES_" + (OP1 + 1) + "_" + (OP2 + 1) + ".txt";
		}


	public void get_FES()
		{
			/* Get the DVPs to use as order parameters in the deltaG free energy calculations */
			if (OP_offset + number_of_points > ROWS) System.err.println("ERROR! The OP offset plus the number of points exceeds the length of the PCs!");
			double[] order_parameter_1 = order_param1.getColumnPackedCopy();
			double[] order_parameter_2 = order_param2.getColumnPackedCopy();
			double[] order_parameter_1_sorted = order_param1.getColumnPackedCopy();
			double[] order_parameter_2_sorted = order_param2.getColumnPackedCopy();

			/* 2D KDE using Gaussian Functions */
			int length = order_parameter_1.length;
			double[] kde_array = new double[2 * length];
			for (int i = 0; i < length; i++)
				{
					kde_array[i + i] = order_parameter_1[i];
					kde_array[i + i + 1] = order_parameter_2[i];
				}
			Arrays.sort(order_parameter_1_sorted);
			Arrays.sort(order_parameter_2_sorted);
			double op1max = order_parameter_1_sorted[length - 1];
			double op2max = order_parameter_2_sorted[length - 1];
			double op1min = order_parameter_1_sorted[0];
			double op2min = order_parameter_2_sorted[0];
			KDE_bounds = new double[] { op1min, op2min, op1max, op2max };
			// System.out.println("The bounds for the KDE are: OP1(" + KDE_bounds[0] + "," + KDE_bounds[2] + "); OP2(" + KDE_bounds[1] + "," + KDE_bounds[3] + ")");

			/* Define the KDE offset and number of points to use in calculating the KDE */
			KDE_offset = 0;
			number_of_points = length;
			if (cellsize == 0)
				{
					double maxDim = Math.max(KDE_bounds[2] - KDE_bounds[0], KDE_bounds[3] - KDE_bounds[1]);
					cellsize = (maxDim / 256);
					// System.out.println("The effective cellsize is: " + cellsize);
				}
			KDE = KernelDensityEstimate2d.compute(kde_array, KDE_offset, number_of_points, KDE_bounds, cellsize, null);
			KDE_bandwidths = DiagonalBandwidthSelector2d.get_Bandwidths();
			// System.out.println("Kernel Bandwidths are: " + nf.format(KDE_bandwidths[0]) + "\t" + nf.format(KDE_bandwidths[3]));
			double[] probabilities = new double[length];
			double[] probabilities_sorted = new double[length];
			for (int i = 0; i < length; i++)
				{
					double prob = KDE.apply(order_parameter_1[i], order_parameter_2[i]);
					probabilities[i] = prob;
					probabilities_sorted[i] = prob;
				}
			Arrays.sort(probabilities_sorted);
			final double prob_max = probabilities_sorted[length - 1];
			final double ln_prob_max = Math.log(prob_max);
			final double KBT = (-0.600); // Units are in kcal/mol, T = 300K (room temp)
			FE = new Matrix(length, 3);
			for (int i = 0; i < length; i++)
				{
					double prob = probabilities[i];
					double ln_prob = Math.log(prob);
					double delta_G = KBT * (ln_prob - ln_prob_max);
					if (delta_G <= 0) delta_G = 0.000; // solves the -0.0000000000000 problem...
					FE.set(i, 0, order_parameter_1[i]);
					FE.set(i, 1, order_parameter_2[i]);
					FE.set(i, 2, delta_G);
				}
			path = file_name_head;
			Matrix_IO.write_Matrix(FE, path, 12, 3);
		}

	public void write_FES_Log()
		{
			LOG = new File(out_directory + "FES_Log_" + description + ".txt");
			try
				{
					LOG_writer = new PrintWriter(new BufferedWriter(new FileWriter(LOG)));
				}
			catch (IOException e)
				{
					System.err.println("Could not create the Job Log file.");
					e.printStackTrace();
				}
			LOG_writer.write("Total Number of Conformations (rows) used to create the DVPs: " + ROWS + "\n");
			LOG_writer.write("Total Number of Modes (cols) used to create the DVPs: " + COLS + "\n");
			LOG_writer.write("Order Parameter 1 (Col# 1): " + (OP1 + 1) + "\n");
			LOG_writer.write("Order Parameter 2 (Col# 2): " + (OP2 + 1) + "\n");
			LOG_writer.write("Number of points to extract from the OPs to use for KDE: " + number_of_points + "\n");
			LOG_writer.write("Offset in selecting points from the order parameters: " + OP_offset + "\n");
			LOG_writer.write("The cellsize used for the 2D grid for KDE: " + nf.format(cellsize) + "\n");
			LOG_writer.write("The bounds for the KDE are: OP1(" + nf.format(KDE_bounds[0]) + "," + nf.format(KDE_bounds[2]) + "); OP2(" + nf.format(KDE_bounds[1]) + "," + nf.format(KDE_bounds[3])
					+ ")" + "\n");
			LOG_writer.write("The kernel bandwidths were: " + nf.format(KDE_bandwidths[0]) + "   " + nf.format(KDE_bandwidths[3]) + "\n\n");
			LOG_writer.write(String.format("%-12s%-12s%-12s%-12s", " ", "OP1", "OP2", "FE"));
			LOG_writer.flush();
			FE.print(LOG_writer, 12, 3);
			LOG_writer.flush();
			LOG_writer.write("\nAnalysis completed: " + DateUtils.now());
			LOG_writer.close();
		}


	public KernelDensityEstimate2d getKDE()
		{
			return KDE;
		}

	public double[] get_KDE_bandwidths()
		{
			return KDE_bandwidths;
		}

	public Matrix get_FE()
		{
			return FE;
		}
}
