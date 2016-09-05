package jed;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import Jama.Matrix;

/**
 * JED class Matrix_Column_Permutation: Randomizes a matrix by permuting its columns.
 * Each column is permuted using the same random 'mapping'.
 * This preserves the orthonormality of eigenvectors.
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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Dr. Charles David
 * 
 */

public class Matrix_Column_Permutation
{

	int ROWS, COLS;
	Matrix data_in, data_permuted;
	File matrix_in, matrix_permuted;
	ArrayList<Integer> indices;
	ArrayList<Double> Components, Components_Permuted;

	/**
	 * Constructor for randomizing a matrix of eigenvectors using a random map permutation applied to each column:
	 * 
	 * @param m
	 *            The matrix whose columns will be permuted.
	 */
	public Matrix_Column_Permutation(Matrix m)
		{

			data_in = m;
			ROWS = data_in.getRowDimension();
			COLS = data_in.getColumnDimension();
			data_permuted = new Matrix(ROWS, COLS);
		}

	/**
	 * Method to establish a random map that is applied to each column of a matrix to permute its elements.
	 * 
	 * @return The randomized matrix
	 */
	public Matrix Get_Random_Permuted_Matrix()
		{

			indices = new ArrayList<Integer>();

			for (int i = 0; i < ROWS; i++)
				{
					indices.add(i);
				}
			Collections.shuffle(indices);

			for (int j = 0; j < COLS; j++)
				{
					int row_index_1 = 0;
					int row_index_2 = ROWS - 1;
					int col_index_1 = j;
					int col_index_2 = j;

					Matrix col = data_in.getMatrix(row_index_1, row_index_2, col_index_1, col_index_2);
					Matrix col_permuted = new Matrix(ROWS, 1);

					for (int i = 0; i < ROWS; i++)
						{
							double element = col.get(i, 0);
							int index = indices.get(i);
							col_permuted.set(index, 0, element);
						}
					data_permuted.setMatrix(row_index_1, row_index_2, col_index_1, col_index_2, col_permuted);
				}
			return data_permuted;
		}
}
