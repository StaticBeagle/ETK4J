package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.MathETK;

public class QRDecomposition {
	protected double[] _data;
	protected final int _rows;
	protected final int _cols;

	protected double[] _rdiag;

	public QRDecomposition(Matrix matrix) {
		final int rows = matrix.getRowCount();
		final int cols = matrix.getColumnCount();
		double[] data = matrix.getArrayCopy();
		_rdiag = new double[cols];

		for (int k = 0; k < cols; ++k) {
			double nrm = 0.0;
			// Compute 2-norm of k-th column without under/overflow.
			for (int i = k; i < rows; ++i) {
				nrm = MathETK.hypot(nrm, data[i * cols + k]);
			}

			if (nrm != 0.0) {
				// Form k-th Householder vector.
				if (data[k * cols + k] < 0) {
					nrm = -nrm;
				}
				for (int i = k; i < rows; i++) {
					data[i * cols + k] /= nrm;
				}
				data[k * cols + k] += 1.0;

				// Apply transformation to remaining columns.
				for (int j = k + 1; j < cols; j++) {
					double s = 0.0;
					for (int i = k; i < rows; i++) {
						s += data[i * cols + k] * data[i * cols + j];
					}
					s = -s / data[k * cols + k];
					for (int i = k; i < rows; i++) {
						data[i * cols + j] += s * data[i * cols + k];
					}
				}
			}
			_rdiag[k] = -nrm;
		}
		_data = data;
		_rows = rows;
		_cols = cols;
	}

	/*
	 * ------------------------ Public Methods ------------------------
	 */

	/**
	 * Is the matrix full rank?
	 * 
	 * @return true if R, and hence A, has full rank.
	 */

	public boolean isFullRank() {
		for (int j = 0; j < _cols; j++) {
			if (_rdiag[j] == 0.0)
				return false;
		}
		return true;
	}

	/**
	 * Return the Householder vectors
	 * 
	 * @return Lower trapezoidal matrix whose columns define the reflections
	 */

	public Matrix getH() {
		Matrix X = new Matrix(_rows, _cols);
		double[] H = X.getArray();
		for (int i = 0; i < _rows; i++) {
			for (int j = 0; j < _cols; j++) {
				if (i >= j) {
					H[i * _cols + j] = _data[i * _cols + j];
				} else {
					H[i * _cols + j] = 0.0;
				}
			}
		}
		return X;
	}

	/**
	 * Return the upper triangular factor
	 * 
	 * @return R
	 */

	public Matrix getR() {
		Matrix X = new Matrix(_cols, _cols);
		double[] R = X.getArray();
		for (int i = 0; i < _cols; i++) {
			for (int j = 0; j < _cols; j++) {
				if (i < j) {
					R[i * _cols + j] = _data[i * _cols + j];
				} else if (i == j) {
					R[i * _cols + j] = _rdiag[i];
				} else {
					R[i * _cols + j] = 0.0;
				}
			}
		}
		return X;
	}
	
	/**
	 * Return the upper triangular factor
	 * 
	 * @return R
	 */

	public Matrix getRT() {
		Matrix X = new Matrix(_cols, _cols);
		double[] R = X.getArray();
		for (int j = 0; j < _cols; ++j) {
			for (int i = 0; i < _cols; ++i) {
				if (i > j) {
					R[i * _cols + j] = _data[j * _cols + i];
				} else if (i == j) {
					R[i * _cols + j] = _rdiag[i];
				} else {
					R[i * _cols + j] = 0.0;
				}
			}
		}
		return X;
	}
	
	/**
	 * Generate and return the (economy-sized) orthogonal factor
	 * 
	 * @return Q
	 */

	public Matrix getQThin() {
		Matrix X = new Matrix(_rows, _cols);
		double[] Q = X.getArray();
		for (int k = _cols - 1; k >= 0; k--) {
			for (int i = 0; i < _rows; i++) {
				Q[i * _cols + k] = 0.0;
			}
			Q[k * _cols + k] = 1.0;
			for (int j = k; j < _cols; j++) {
				if (_data[k * _cols + k] != 0) {
					double s = 0.0;
					for (int i = k; i < _rows; i++) {
						s += _data[i * _cols + k] * Q[i * _cols + j];
					}
					s = -s / _data[k * _cols + k];
					for (int i = k; i < _rows; i++) {
						Q[i * _cols + j] += s * _data[i * _cols + k];
					}
				}
			}
		}
		return X;
	}

	/**
	 * Generate and return the unitary orthogonal factor
	 * 
	 * @return Q
	 */
	public Matrix getQ() {
		// Compute Q = Q * I
		Matrix Q = Matrix.identity(_rows, _rows);
		double[] X = Q.getArray();
		int mr = Math.min(_rows, _cols);
		for (int k = mr - 1; k >= 0; --k) {
			for (int j = _rows - 1; j >= 0; --j) {
				double s = 0.0;
				for (int i = k; i < _rows; i++) {
					s += _data[i * _cols + k] * X[i * _rows + j];
				}
				s = -s / _data[k * _cols + k];
				for (int i = k; i < _rows; i++) {
					X[i * _rows + j] += s * _data[i * _cols + k];
				}
			}
		}
		return Q;
	}
	
	public Matrix QmultiplyX(Matrix X) {
		int nx = X.getColumnCount();
		double[] Q = X.getArray();
		int mr = Math.min(_rows, _cols);
		for (int k = mr - 1; k >= 0; --k) {
			for (int j = nx - 1; j >= 0; --j) {
				double s = 0.0;
				for (int i = k; i < _rows; i++) {
					s += _data[i * _cols + k] * Q[i * nx + j];
				}
				s = -s / _data[k * _cols + k];
				for (int i = k; i < _rows; i++) {
					Q[i * nx + j] += s * _data[i * _cols + k];
				}
			}
		}
		return new Matrix(Q, X.getRowCount(), X.getColumnCount());
	}
	
	/**
	 * Generate and return the transpose of the orthogonal factor
	 * 
	 * @return transpose(Q)
	 */
	public Matrix getQT() {
		// Compute Q = Q * I
		Matrix Q = Matrix.identity(_rows, _rows);
		double[] X = Q.getArray();
		int mr = Math.min(_rows, _cols);
		for (int k = 0; k < mr; ++k) {
			for (int j = 0; j < _rows; ++j) {
				double s = 0.0;
				for (int i = k; i < _rows; i++) {
					s += _data[i * _cols + k] * X[i * _rows + j];
				}
				s = -s / _data[k * _cols + k];
				for (int i = k; i < _rows; i++) {
					X[i * _rows + j] += s * _data[i * _cols + k];
				}
			}
		}
		return Q;
	}

	/**
	 * Least squares solution of A*X = B
	 * 
	 * @param B
	 *            A Matrix with as many rows as A and any number of columns.
	 * @return X that minimizes the two norm of Q*R*X-B.
	 * @exception IllegalArgumentException
	 *                Matrix row dimensions must agree.
	 * @exception RuntimeException
	 *                Matrix is rank deficient.
	 */

	public Matrix solve(Matrix B) {
		if (B.getRowCount() != _rows) {
			throw new IllegalArgumentException("Matrix row dimensions must agree.");
		}
		if (!this.isFullRank()) {
			throw new RuntimeException("Matrix is rank deficient.");
		}

		// Copy right hand side
		int nx = B.getColumnCount();
		double[] X = B.getArrayCopy();

		// Compute Y = transpose(Q)*B
		for (int k = 0; k < _cols; k++) {
			for (int j = 0; j < nx; j++) {
				double s = 0.0;
				for (int i = k; i < _rows; i++) {
					s += _data[i * _cols + k] * X[i * nx + j];
				}
				s = -s / _data[k * _cols + k];
				for (int i = k; i < _rows; i++) {
					X[i * nx + j] += s * _data[i * _cols + k];
				}
			}
		}
		// Solve R*X = Y;
		for (int k = _cols - 1; k >= 0; k--) {
			for (int j = 0; j < nx; j++) {
				X[k * nx + j] /= _rdiag[k];
			}
			for (int i = 0; i < k; i++) {
				for (int j = 0; j < nx; j++) {
					X[i * nx + j] -= X[k * nx + j] * _data[i * _cols + k];
				}
			}
		}
		return (new Matrix(X, _cols, nx).subMatrix(0, _cols - 1, 0, nx - 1));
	}
	
//	public Matrix solveTranspose(Matrix B) {
////		if (B.getRowCount() != _rows) {
////			throw new IllegalArgumentException("Matrix row dimensions must agree.");
////		}
//		if (!this.isFullRank()) {
//			throw new RuntimeException("Matrix is rank deficient.");
//		}
//
//		// Copy right hand side
//		int nx = B.getColumnCount();
//		double[] X = B.getArrayCopy();
//		
//		// Solve RT*X = Y;
//		for (int k = _cols - 1; k >= 0; k--) {
//			for (int j = 0; j < nx; j++) {
//				X[k * nx + j] /= _rdiag[k];
//			}
//			for (int i = 0; i < k; i++) {
//				for (int j = 0; j < nx; j++) {
//					X[i * nx + j] -= X[k * nx + j] * _data[i * _cols + k];
//				}
//			}
//		}
//
//
//		int mr = Math.min(_rows, _cols);
//		for (int k = mr - 1; k >= 0; --k) {
//			for (int j = nx - 1; j >= 0; --j) {
//				double s = 0.0;
//				for (int i = k; i < _rows; i++) {
//					s += _data[i * _cols + k] * X[i * nx + j];
//				}
//				s = -s / _data[k * _cols + k];
//				for (int i = k; i < _rows; i++) {
//					X[i * nx + j] += s * _data[i * _cols + k];
//				}
//			}
//		}
//		return new Matrix(X, B.getRowCount(), B.getColumnCount());
//
//		//return (new Matrix(X, _cols, nx).subMatrix(0, _cols - 1, 0, nx - 1));
//	}

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
