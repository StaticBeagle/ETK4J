package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

final class SolverUtils {
	private SolverUtils() {}

	public static class LU {
		protected double[][] _data;
		protected int _cols;

		protected int _pivotsign;
		protected int[] _pivot;

		public LU(double[][] matrix) {
			int rows = matrix.length;
			int cols = _cols = matrix[0].length;
			_data = matrix;

			_pivot = new int[rows];
			for (int i = 0; i < rows; i++) {
				_pivot[i] = i;
			}
			_pivotsign = 1;

			double[] LUrow;
			double[] LUcol = new double[rows];

			// Begin the outer loop
			for (int j = 0; j < cols; j++) {
				// Copy the j-th column to localize references.
				for (int i = 0; i < rows; i++) {
					LUcol[i] = _data[i][j];
				}
				// Apply previous transformations
				for (int i = 0; i < rows; i++) {
					LUrow = _data[i];

					int maxel = Math.min(i, j);
					double s = 0.0;
					for (int k = 0; k < maxel; k++) {
						s += LUrow[k] * LUcol[k];
					}

					LUrow[j] = LUcol[i] -= s;
				}

				// Find pivot and swap if needed
				int p = j;
				for (int i = j + 1; i < rows; i++) {
					if (Double.compare(Math.abs(LUcol[i]), Math.abs(LUcol[p])) > 0) {
						p = i;
					}
				}
				if (p != j) {
					for (int k = 0; k < cols; k++) {
						double temp = _data[p][k];
						_data[p][k] = _data[j][k];
						_data[j][k] = temp;
					}
					int temp = _pivot[p];
					_pivot[p] = _pivot[j];
					_pivot[j] = temp;
				}
				// Wrapping up
				if (j < rows & Double.compare(_data[j][j], 0.0) != 0) {
					for (int i = j + 1; i < rows; i++) {
						_data[i][j] /= _data[j][j];
					}
				}
			}
		}

		public double det() {
			double det = _pivotsign;
			for (int i = 0; i < _cols; i++) {
				det *= _data[i][i];
			}
			return det;
		}

		public void solve(final double[] vector, double[] solution) {
			int rows = vector.length;
			// Shuffle rows to match pivot vector
			for (int i = 0; i < rows; i++) {
				solution[i] = vector[_pivot[i]];
			}

			// Solve Ly = b
			// Forward substitution
			for (int k = 0; k < rows; k++) {
				for (int i = k + 1; i < rows; i++) {
					solution[i] -= solution[k] * _data[i][k];
				}
			}

			// Solve Ux = y
			// Back substitution
			for (int k = rows - 1; k >= 0; k--) {
				solution[k] /= _data[k][k];
				for (int i = 0; i < k; i++) {
					solution[i] -= solution[k] * _data[i][k];
				}
			}
		}
	}
}
