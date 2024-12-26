package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.util.DoubleArrays;

/**
 * Cholesky Decomposition.
 * <P>
 * For a symmetric, positive definite matrix A, the Cholesky decomposition is an
 * lower triangular matrix L so that A = L*L'.
 * <P>
 * If the matrix is not symmetric or positive definite, the constructor returns
 * a partial decomposition and sets an internal flag that may be queried by the
 * isSPD() method.
 */

public class CholeskyDecompositionDense extends CholeskyDecomposition<MatrixDense>{

	/*
	 * ------------------------ Class variables ------------------------
	 */

	/**
	 * Array for internal storage of decomposition.
	 * 
	 * @serial internal array storage.
	 */
	private double[][] L;

	/**
	 * Row and column dimension (square matrix).
	 * 
	 * @serial matrix dimension.
	 */
	private int _n;

	/**
	 * Symmetric and positive definite flag.
	 * 
	 * @serial is symmetric and positive definite flag.
	 */
	private boolean isspd;

	/*
	 * ------------------------ Constructor ------------------------
	 */

	/**
	 * Cholesky algorithm for symmetric and positive definite matrix. Structure
	 * to access L and isspd flag.
	 * 
	 * @param Arg
	 *            Square, symmetric matrix.
	 */

	public CholeskyDecompositionDense(MatrixDense Arg) {
		super(Arg);
		// Initialize.
		double[] A = Arg.getArray();
		final int m = Arg.getRowCount();
		final int n = Arg.getColumnCount();
		_n = m;
		L = new double[m][m];
		isspd = (n == m);
		// Main loop.
		for (int j = 0; j < m; j++) {
			double[] Lrowj = L[j];
			double d = 0.0;
			for (int k = 0; k < j; k++) {
				double[] Lrowk = L[k];
				double s = 0.0;
				for (int i = 0; i < k; i++) {
					s += Lrowk[i] * Lrowj[i];
				}
				Lrowj[k] = s = (A[j * n + k] - s) / L[k][k];
				d = d + s * s;
				isspd = isspd & (A[k * n + j] == A[j * n + k]);
			}
			d = A[j * n + j] - d;
			isspd = isspd & (d > 0.0);
			L[j][j] = Math.sqrt(Math.max(d, 0.0));
			for (int k = j + 1; k < m; k++) {
				L[j][k] = 0.0;
			}
		}
	}
	// TODO change from Lower to upper R triangular matrix
	/*
	 * ------------------------ Temporary, experimental code.
	 * ------------------------ *\
	 * 
	 * \** Right Triangular Cholesky Decomposition. <P> For a symmetric,
	 * positive definite matrix A, the Right Cholesky decomposition is an upper
	 * triangular matrix R so that A = R'*R. This constructor computes R with
	 * the Fortran inspired column oriented algorithm used in LINPACK and
	 * MATLAB. In Java, we suspect a row oriented, lower triangular
	 * decomposition is faster. We have temporarily included this constructor
	 * here until timing experiments confirm this suspicion. \
	 * 
	 * \** Array for internal storage of right triangular decomposition. **\
	 * private transient double[][] R;
	 * 
	 * \** Cholesky algorithm for symmetric and positive definite matrix.
	 * 
	 * @param A Square, symmetric matrix.
	 * 
	 * @param rightflag Actual value ignored.
	 * 
	 * @return Structure to access R and isspd flag. \
	 * 
	 * public CholeskyDecomposition (Matrix Arg, int rightflag) { // Initialize.
	 * double[][] A = Arg.getArray(); n = Arg.getColumnDimension(); R = new
	 * double[n][n]; isspd = (Arg.getColumnDimension() == n); // Main loop. for
	 * (int j = 0; j < n; j++) { double d = 0.0; for (int k = 0; k < j; k++) {
	 * double s = A[k][j]; for (int i = 0; i < k; i++) { s = s -
	 * R[i][k]*R[i][j]; } R[k][j] = s = s/R[k][k]; d = d + s*s; isspd = isspd &
	 * (A[k][j] == A[j][k]); } d = A[j][j] - d; isspd = isspd & (d > 0.0);
	 * R[j][j] = Math.sqrt(Math.max(d,0.0)); for (int k = j+1; k < n; k++) {
	 * R[k][j] = 0.0; } } }
	 * 
	 * \** Return upper triangular factor.
	 * 
	 * @return R \
	 * 
	 * public Matrix getR () { return new Matrix(R,n,n); }
	 * 
	 * \* ------------------------ End of temporary code.
	 * ------------------------
	 */

	/*
	 * ------------------------ Public Methods ------------------------
	 */

	/**
	 * Is the matrix symmetric and positive definite?
	 * 
	 * @return true if A is symmetric and positive definite.
	 */

	public boolean isSPD() {
		return isspd;
	}

	/**
	 * Return triangular factor.
	 * 
	 * @return L
	 */

	public MatrixDense getL() {
		return new MatrixDense(DoubleArrays.flatten(L), _n, _n);
	}

	/**
	 * Solve A*X = B
	 * 
	 * @param B
	 *            A Matrix with as many rows as A and any number of columns.
	 * @return X so that L*L'*X = B
	 * @exception IllegalArgumentException
	 *                Matrix row dimensions must agree.
	 * @exception RuntimeException
	 *                Matrix is not symmetric positive definite.
	 */

	public MatrixDense solve(MatrixDense B) {
		final int n = _n;
		if (B.getRowCount() != n) {
			throw new IllegalArgumentException("Matrix row dimensions must agree.");
		}
		if (!isspd) {
			throw new RuntimeException("Matrix is not symmetric positive definite.");
		}

		// Copy right hand side.
		double[] X = B.getArrayCopy();
		int nx = B.getColumnCount();

		// Solve L*Y = B;
		for (int k = 0; k < n; k++) {
			for (int j = 0; j < nx; j++) {
				for (int i = 0; i < k; i++) {
					X[k * nx + j] -= X[i * nx + j] * L[k][i];
				}
				X[k * nx + j] /= L[k][k];
			}
		}

		// Solve L'*X = Y;
		for (int k = n - 1; k >= 0; k--) {
			for (int j = 0; j < nx; j++) {
				for (int i = k + 1; i < n; i++) {
					X[k * nx + j] -= X[i * nx + j] * L[i][k];
				}
				X[k * nx + j] /= L[k][k];
			}
		}

		return new MatrixDense(X, n, nx);
	}
}
