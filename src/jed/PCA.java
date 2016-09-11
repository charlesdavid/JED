package jed;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 * JED class PCA: Constructs the COV(Q), CORR(R), and PCORR(P) matrices. Provides a hook for the Jama eigenvalue decomposition. Note: The assumption is that ROWS are variables and
 * COLS are observations Copyright (C) 2012 Dr. Charles David
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
 */

public class PCA
{
	int ROWS, COLS;
	Matrix data_matrix, data_means, data_sigmas;

	/* ****************************** CONSTRUCTOR ********************************************************************************* */

	/**
	 * Constructor to perform PCA on data: Note: The data matrix takes ROWS as variables and COLS as instances.
	 * 
	 * @param data
	 *            The data matrix A
	 */
	public PCA(Matrix data)
		{
			data_matrix = data;
			COLS = data.getColumnDimension();
			ROWS = data.getRowDimension();

		}

	/* ****************************** METHODS ********************************************************************************* */

	/**
	 * Method to calculate the covariance matrix: AA(T) Note: The data are ROW CENTERED before the calculation.
	 * 
	 * @return Returns the Covariance Matrix, Q
	 */
	Matrix get_covariance_matrix()
		{

			Row_Center_Data rcd = new Row_Center_Data(data_matrix);
			data_matrix = null;
			System.gc();

			Matrix X_vectors = rcd.get_row_centered_data();
			data_means = rcd.get_variable_means();
			data_sigmas = rcd.get_variable_sigmas();

			Matrix Q = X_vectors.times(X_vectors.transpose());
			Q.timesEquals(Math.pow((COLS - 1), -1));

			X_vectors = null;
			System.gc();
			return Q;
		}

	/**
	 * Method to calculate the covariance matrix (without matrix multiplication): Note: The data are ROW CENTERED before the calculation.
	 * 
	 * @return Returns the Covariance Matrix, Q
	 */
	Matrix get_covariance_matrix_elegant()
		{

			Row_Center_Data rcd = new Row_Center_Data(data_matrix);
			data_means = rcd.get_variable_means();
			data_sigmas = rcd.get_variable_sigmas();

			Matrix Q = new Matrix(ROWS, ROWS);

			for (int i = 0; i < ROWS; i++)
			{
				double s = data_sigmas.get(i, 0);
				double var = (s * s);
				Q.set(i, i, var);
				double mean_X = data_means.get(i, 0);
				Matrix var_X = data_matrix.getMatrix(i, i, 0, COLS - 1);
				for (int j = i + 1; j < ROWS; j++)
				{
					Matrix var_Y = data_matrix.getMatrix(j, j, 0, COLS - 1);
					double[] var_XY = (var_X.arrayTimes(var_Y)).getRowPackedCopy();
					double mean_Y = data_means.get(j, 0);
					double mean_XY = Descriptive_Stats.get_mean(var_XY);
					double cov = (mean_XY - (mean_X * mean_Y));
					Q.set(i, j, cov);
					Q.set(j, i, cov);
				}
			}
			data_matrix = null;
			System.gc();
			return Q;
		}

	/**
	 * Method to calculate the covariance matrix: AA(T) Note: The data are NOT CENTERED before the calculation.
	 * 
	 * @return Returns the Covariance Matrix, Q
	 */
	Matrix get_covariance_matrix_NO_CENTERING()
		{

			Matrix Q = data_matrix.times(data_matrix.transpose());
			data_matrix = null;
			System.gc();
			Q.timesEquals(Math.pow((COLS - 1), -1));
			return Q;
		}

	/**
	 * Method to calculate the correlation matrix from a covariance matrix.
	 * 
	 * @return Returns the Correlation Matrix, R
	 */
	Matrix get_R_from_Q(Matrix Cov_Matrix)
		{

			int num_Vars = Cov_Matrix.getColumnDimension();
			Matrix R = new Matrix(num_Vars, num_Vars);
			for (int i = 0; i < num_Vars; i++)
			{
				double sigma = Math.sqrt(Cov_Matrix.get(i, i));
				R.set(i, i, 1d);
				for (int j = 0; j < i; j++)
				{
					double element = Cov_Matrix.get(i, j) / (sigma * Math.sqrt(Cov_Matrix.get(j, j)));
					R.set(j, i, element);
					R.set(i, j, element);
				}
			}
			return R;
		}

	/**
	 * Method to calculate the correlation matrix from the data. Note: The data are ROW CENTERED before the calculation.
	 * 
	 * @return Returns the Correlation Matrix, R
	 */
	Matrix get_correlation_matrix()
		{
			Matrix COV = get_covariance_matrix();
			Matrix R = get_R_from_Q(COV);
			return R;
		}

	/**
	 * Method to calculate the correlation matrix from the data. Note: The data are NOT CENTERED before the calculation.
	 * 
	 * @return Returns the Correlation Matrix, R
	 */
	Matrix get_correlation_matrix_NO_CENTERING()
		{
			Matrix COV = get_covariance_matrix_NO_CENTERING();
			Matrix R = get_R_from_Q(COV);
			return R;
		}

	/**
	 * Method to calculate the partial-correlation matrix from a precision matrix.
	 * 
	 * @return Returns the Partial Correlation Matrix, P_CORR
	 */
	static Matrix get_partial_correlation_matrix(Matrix precision)
		{
			int num_Vars = precision.getColumnDimension();
			Matrix P_CORR = new Matrix(num_Vars, num_Vars);
			for (int i = 0; i < num_Vars; i++)
			{
				double x = precision.get(i, i);
				P_CORR.set(i, i, -1d);
				for (int j = 0; j < i; j++)
				{
					double y = precision.get(j, j);
					double r = precision.get(i, j);
					double p = (-r / (Math.sqrt(x * y)));
					P_CORR.set(j, i, p);
					P_CORR.set(i, j, p);
				}
			}

			return P_CORR;
		}

	/**
	 * Method to calculate the eigenvalue decomposition of a matrix. Note: The matrix argument should be a positive-definite matrix for real eigenvalues. Based in the JAMA library
	 * 
	 * @param q
	 *            The matrix to factor
	 * @return The eigenvalue decomposition holding the eigenvalues and eigenvectors.
	 */
	static EigenvalueDecomposition get_eigenvalue_decomposition(Matrix q)
		{

			EigenvalueDecomposition evd = new EigenvalueDecomposition(q);

			return evd;
		}

	/**
	 * @return Returns the means of the centered variables.
	 */
	public Matrix getData_means()
		{
			return data_means;
		}

	/**
	 * @return Returns the standard deviations of the centered variables.
	 */
	public Matrix getData_sigmas()
		{
			return data_sigmas;
		}
}