package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

import java.util.Arrays;

public class ComplexLUDecomposition {
	protected Complex[] _data;
	protected final int _rows;
	protected final int _cols;

	protected int _pivotsign = 1;
	protected int[] _pivot;

	public ComplexLUDecomposition(ComplexMatrix matrix) {
		final int rows = matrix.getRowCount();
		final int cols = matrix.getColumnCount();
		Complex[] data = matrix.getArrayCopy();

		_pivot = new int[rows];
		for (int i = 0; i < rows; i++) {
			_pivot[i] = i;
		}

		Complex[] LUcol = new Complex[rows];

		// Begin the outer loop
		for (int j = 0; j < cols; j++) {
			// Copy the j-th column to localize references.
			for (int i = 0; i < rows; i++) {
				LUcol[i] = data[i * cols + j];
			}
			// Apply previous transformations
			for (int i = 0; i < rows; i++) {
				int maxel = Math.min(i, j);
				Complex s = new Complex();
				for (int k = 0; k < maxel; k++) {
					s.addEquals(data[i * cols + k].multiply(LUcol[k]));
				}
				LUcol[i].subtractEquals(s);
				data[i * cols + j] = LUcol[i];
			}

			// Find pivot and swap if needed
			int p = j;
			for (int i = j + 1; i < rows; i++) {
				if (LUcol[i].abs() > LUcol[p].abs()) {
					p = i;
				}
			}
			if (p != j) {
				for (int k = 0; k < cols; k++) {
					Complex temp = data[p * cols + k];
					data[p * cols + k] = data[j * cols + k];
					data[j * cols + k] = temp;
				}
				int temp = _pivot[p];
				_pivot[p] = _pivot[j];
				_pivot[j] = temp;
				_pivotsign = -_pivotsign;
			}
			// Wrapping up
			if (j < rows & !data[j * cols + j].equals(new Complex())) {
				for (int i = j + 1; i < rows; i++) {
					data[i * cols + j].divideEquals(data[j * cols + j]);
				}
			}
		}
		_data = data;
		_rows = rows;
		_cols = cols;
	}

	public boolean isNonSingular() {
		for (int j = 0; j < _cols; ++j) {
			if (_data[j * _cols + j].equals(new Complex())) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Return lower triangular factor
	 * 
	 * @return L
	 */

	public ComplexMatrix getL() {
		final int rows = _rows;
		final int cols = _cols;
		Complex[] L = new Complex[rows * cols];
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				if (i > j) {
					L[i * cols + j] = _data[i * cols + j];
				} else if (i == j) {
					L[i * cols + j] = Complex.fromReal(1.0);
				} else {
					L[i * cols + j] = new Complex();
				}
			}
		}
		return new ComplexMatrix(L, rows, cols);
	}

	/**
	 * Return upper triangular factor
	 * 
	 * @return U
	 */

	public ComplexMatrix getU() {
		final int rows = _rows;
		final int cols = _cols;
		Complex[] U = new Complex[rows * cols];
		for (int i = 0; i < cols; i++) {
			for (int j = 0; j < cols; j++) {
				if (i <= j) {
					U[i * cols + j] = _data[i * cols + j];
				} else {
					U[i * cols + j] = new Complex();
				}
			}
		}
		return new ComplexMatrix(U, rows, cols);
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
		final int rows = _rows;
		double[] vals = new double[rows];
		for (int i = 0; i < rows; i++) {
			vals[i] = (double) _pivot[i];
		}
		return vals;
	}

	public Complex det() {
		if (_rows != _cols) {
			throw new IllegalArgumentException("Matrix must be square.");
		}
		Complex det = Complex.fromReal(_pivotsign);
		for (int i = 0; i < _cols; i++) {
			det.multiplyEquals(this._data[i * _cols + i]);
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

	public ComplexMatrix solve(Matrix B) {
		if (B.getRowCount() != _rows) {
			throw new IllegalArgumentException("Matrix row dimensions must agree.");
		}
		if (!this.isNonSingular()) {
			throw new RuntimeException("Matrix is singular.");
		}

		// Copy right hand side with pivoting
		int nx = B.getColumnCount();
		ComplexMatrix Xmat = new ComplexMatrix(B.subMatrix(_pivot, 0, nx - 1));
		Complex[] X = Xmat.getArray();

		final int cols = _cols;

		// Solve L * Y = B(_pivot,:)
		for (int k = 0; k < cols; ++k) {
			for (int i = k + 1; i < cols; ++i) {
				for (int j = 0; j < nx; ++j) {
					X[i * nx + j].subtractEquals(X[k * nx + j].multiply(_data[i * cols + k]));
				}
			}
		}
		// Solve U * X = Y;
		for (int k = cols - 1; k >= 0; --k) {
			for (int j = 0; j < nx; ++j) {
				X[k * nx + j].divideEquals(_data[k * cols + k]);
			}
			for (int i = 0; i < k; ++i) {
				for (int j = 0; j < nx; ++j) {
					X[i * nx + j].subtractEquals(X[k * nx + j].multiply(_data[i * cols + k]));
				}
			}
		}
		return Xmat;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < _rows * _cols; ++i) {
			if (i > 0 && i % _cols == 0) {
				sb.append(System.lineSeparator());
			}
			sb.append(String.format("%.4f", _data[i])).append(" ");
		}
		return sb.toString();
	}
}