package com.wildbitsfoundry.etk4j.math.linearalgebra;

import java.util.Arrays;

public class LUDecompositionDense extends LUDecomposition<MatrixDense> {
	protected double[] _data;

	protected int _pivotsign = 1;
	protected int[] _pivot;

	private MatrixDense L;
	private MatrixDense U;

	//private

	public LUDecompositionDense(MatrixDense matrix) {
		super(matrix);
		final int rows = matrix.getRowCount();
		final int cols = matrix.getColumnCount();
		double[] data = matrix.getArrayCopy();

		_pivot = new int[rows];
		for (int i = 0; i < rows; i++) {
			_pivot[i] = i;
		}

		double[] LUcol = new double[rows];

		// Begin the outer loop
		for (int j = 0; j < cols; j++) {
			// Copy the j-th column to localize references.
			for (int i = 0; i < rows; i++) {
				LUcol[i] = data[i * cols + j];
			}
			// Apply previous transformations
			for (int i = 0; i < rows; i++) {
				int maxel = Math.min(i, j);
				double s = 0.0;
				for (int k = 0; k < maxel; k++) {
					s += data[i * cols + k] * LUcol[k];
				}

				data[i * cols + j] = LUcol[i] -= s;
			}

			// Find pivot and swap if needed
			int p = j;
			for (int i = j + 1; i < rows; i++) {
				if (Math.abs(LUcol[i]) > Math.abs(LUcol[p])) {
					p = i;
				}
			}
			if (p != j) {
				for (int k = 0; k < cols; k++) {
					double temp = data[p * cols + k];
					data[p * cols + k] = data[j * cols + k];
					data[j * cols + k] = temp;
				}
				int temp = _pivot[p];
				_pivot[p] = _pivot[j];
				_pivot[j] = temp;
				_pivotsign = -_pivotsign;
			}
			// Wrapping up
			if (j < rows & data[j * cols + j] != 0.0) {
				for (int i = j + 1; i < rows; i++) {
					data[i * cols + j] /= data[j * cols + j];
				}
			}
		}
		_data = data;
	}

	public boolean isSingular() {
		for (int j = 0; j < cols; ++j) {
			if (_data[j * cols + j] == 0) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Return lower triangular factor
	 * 
	 * @return L
	 */

	public MatrixDense getL() {
		if(this.L != null) {
			return this.L;
		}
		final int rows = this.rows;
		final int cols = this.cols;
		double[] L = new double[rows * cols];
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				if (i > j) {
					L[i * cols + j] = _data[i * cols + j];
				} else if (i == j) {
					L[i * cols + j] = 1.0;
				} else {
					L[i * cols + j] = 0.0;
				}
			}
		}
		this.L = new MatrixDense(L, rows, cols);
		return this.L;
	}

	/**
	 * Return upper triangular factor
	 * 
	 * @return U
	 */

	public MatrixDense getU() {
		if(this.U != null) {
			return this.U;
		}
		final int rows = this.rows;
		final int cols = this.cols;
		double[] U = new double[rows * cols];
		for (int i = 0; i < cols; i++) {
			for (int j = 0; j < cols; j++) {
				if (i <= j) {
					U[i * cols + j] = _data[i * cols + j];
				} else {
					U[i * cols + j] = 0.0;
				}
			}
		}
		this.U = new MatrixDense(U, rows, cols);
		return this.U;
	}

	/**
	 * Return pivot permutation vector
	 * 
	 * @return piv
	 */

	public int[] getPivot() {
		return Arrays.copyOf(_pivot, _pivot.length);
	}

	/**
	 * Return pivot permutation vector as a one-dimensional double array
	 * 
	 * @return (double) piv
	 */

	public double[] getPivotAsDouble() {
		final int rows = this.rows;
		double[] vals = new double[rows];
		for (int i = 0; i < rows; i++) {
			vals[i] = _pivot[i];
		}
		return vals;
	}

	public double det() {
		if (this.rows != this.cols) {
			throw new IllegalArgumentException("Matrix must be square.");
		}
		double det = _pivotsign;
		for (int i = 0; i < this.cols; i++) {
			det *= this._data[i * this.cols + i];
		}
		return det;
	}

	/**
	 * Solve A*X = B
	 * 
	 * @param B
	 *            A Matrix with as many rows as A and any number of columns.
	 * @return X so that L*U*X = B(piv,:)
	 * @exception IllegalArgumentException
	 *                Matrix row dimensions must agree.
	 * @exception RuntimeException
	 *                Matrix is singular.
	 */

	public MatrixDense solve(MatrixDense B) {
		if (B.getRowCount() != this.rows) {
			throw new IllegalArgumentException("Matrix row dimensions must agree.");
		}
		if (this.isSingular()) {
			throw new RuntimeException("Matrix is singular.");
		}

		// Copy right-hand side with pivoting
		int nx = B.getColumnCount();
		MatrixDense Xmat = B.subMatrix(_pivot, 0, nx - 1);
		double[] X = Xmat.getArray();

		final int cols = this.cols;

		// Solve L * Y = B(_pivot,:)
		for (int k = 0; k < cols; ++k) {
			for (int i = k + 1; i < cols; ++i) {
				for (int j = 0; j < nx; ++j) {
					X[i * nx + j] -= X[k * nx + j] * _data[i * cols + k];
				}
			}
		}
		// Solve U * X = Y;
		for (int k = cols - 1; k >= 0; --k) {
			for (int j = 0; j < nx; ++j) {
				X[k * nx + j] /= _data[k * cols + k];
			}
			for (int i = 0; i < k; ++i) {
				for (int j = 0; j < nx; ++j) {
					X[i * nx + j] -= X[k * nx + j] * _data[i * cols + k];
				}
			}
		}
		return Xmat;
	}

	// TODO write test
	@Override
	public MatrixDense solve(double[] b) {
		return solve(new MatrixDense(b, b.length));
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < this.rows * this.cols; ++i) {
			if (i > 0 && i % this.cols == 0) {
				sb.append(System.lineSeparator());
			}
			sb.append(String.format("%.4f", _data[i])).append(" ");
		}
		return sb.toString();
	}
}