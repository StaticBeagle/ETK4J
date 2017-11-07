package com.wildbitsfoundry.etk4j.math.linearalgebra;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.constants.ETKConstants;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.util.NumArrays;

public class Matrix {
	private double[] _data;
	private int _rows;
	private int _cols;

	public Matrix(int rows, int cols) {
		_rows = rows;
		_cols = cols;

		this._data = new double[rows * cols];
	}

	public Matrix(double[] data, int rows) {
		_rows = rows;
		_cols = data.length / rows;
		int dim = _rows * _cols;
		_data = new double[dim];
		System.arraycopy(data, 0, _data, 0, dim);

	}

	public Matrix(double[][] data) {
		_rows = data.length;
		_cols = data[0].length;
		_data = NumArrays.flatten(data);
	}

	public Matrix(double[][] data, int rows, int cols) {
		_rows = rows;
		_cols = cols;
		_data = NumArrays.flatten(data);
	}

	public Matrix(double[] data, int rows, int cols) {
		_rows = rows;
		_cols = cols;
		_data = data;
	}

	public Matrix(Matrix matrix) {
		_rows = matrix._rows;
		_cols = matrix._cols;
		_data = new double[_rows * _cols];
		System.arraycopy(matrix._data, 0, this._data, 0, this._rows * this._cols);
	}

	public Matrix(int rows, int cols, double val) {
		_rows = rows;
		_cols = cols;
		_data = new double[_rows * _cols];
		Arrays.fill(_data, val);
	}

	/***
	 * Get submatrix
	 * 
	 * @param row0
	 *            Initial row index
	 * @param row1
	 *            Final row index
	 * @param col0
	 *            Initial column index
	 * @param col1
	 *            Final column index
	 * @return A(row0 : row1, col0 : col1)
	 */
	public Matrix subMatrix(int row0, int row1, int col0, int col1) {
		int rowDim = row1 - row0 + 1;
		int colDim = col1 - col0 + 1;
		double[] data = new double[rowDim * colDim];
		for (int i = 0; i < rowDim; ++i) {
			for (int j = 0; j < colDim; ++j) {
				data[i * colDim + j] = _data[(i + row0) * _cols + (j + col0)];
			}
		}
		return new Matrix(data, rowDim, colDim);
	}

	/***
	 * Get submatrix
	 * 
	 * @param rows
	 *            Array of row indices
	 * @param col0
	 *            Initial column index
	 * @param col1
	 *            Final column index
	 * @return A(rows(:), col0 : col1)
	 */
	public Matrix subMatrix(int[] rows, int col0, int col1) {
		int rowDim = rows.length;
		int colDim = col1 - col0 + 1;
		double[] data = new double[rowDim * colDim];
		for (int i = 0; i < rowDim; ++i) {
			for (int j = 0; j < colDim; ++j) {
				data[i * colDim + j] = _data[rows[i] * _cols + (j + col0)];
			}
		}
		return new Matrix(data, rowDim, colDim);
	}

	/***
	 * Get submatrix
	 * 
	 * @param row0
	 *            Initial row index
	 * @param row1
	 *            Final row index
	 * @param cols
	 *            Array of column indices
	 * @return A(row0 : row1, cols(:))
	 */
	public Matrix subMatrix(int row0, int row1, int[] cols) {
		int rowDim = row1 - row0 + 1;
		int colDim = cols.length;
		;
		double[] data = new double[rowDim * colDim];
		for (int i = 0; i < rowDim; ++i) {
			for (int j = 0; j < colDim; ++j) {
				data[i * colDim + j] = _data[(i + row0) * _cols + cols[j]];
			}
		}
		return new Matrix(data, rowDim, colDim);
	}

	/***
	 * Get submatrix
	 * 
	 * @param rows
	 *            Array or row indices
	 * @param cols
	 *            Array of column indices
	 * @return A(rows(:), cols(:))
	 */
	public Matrix subMatrix(int[] rows, int[] cols) {
		int rowDim = rows.length;
		int colDim = cols.length;
		double[] data = new double[rowDim * colDim];
		for (int i = 0; i < rowDim; ++i) {
			for (int j = 0; j < colDim; ++j) {
				data[i * colDim + j] = _data[rows[i] * _cols + cols[j]];
			}
		}
		return new Matrix(data, rowDim, colDim);
	}

	public double get(int i, int j) {
		return _data[i * _cols + j];
	}

	public void set(int i, int j, double val) {
		_data[i * _cols + j] = val;
	}

	public double det() {
		return new LUDecomposition(this).det();
	}

	/**
	 * Set a submatrix.
	 * 
	 * @param i0
	 *            Initial row index
	 * @param i1
	 *            Final row index
	 * @param j0
	 *            Initial column index
	 * @param j1
	 *            Final column index
	 * @param X
	 *            A(i0:i1,j0:j1)
	 * @exception ArrayIndexOutOfBoundsException
	 *                Submatrix indices
	 */

	public void setMatrix(int i0, int i1, int j0, int j1, Matrix X) {
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = j0; j <= j1; j++) {
					_data[i * _cols + j] = X.get(i - i0, j - j0);
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
	}

	/**
	 * Set a submatrix.
	 * 
	 * @param r
	 *            Array of row indices.
	 * @param c
	 *            Array of column indices.
	 * @param X
	 *            A(r(:),c(:))
	 * @exception ArrayIndexOutOfBoundsException
	 *                Submatrix indices
	 */

	public void setMatrix(int[] r, int[] c, Matrix X) {
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = 0; j < c.length; j++) {
					_data[r[i] * _cols + c[j]] = X.get(i, j);
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
	}

	/**
	 * Set a submatrix.
	 * 
	 * @param r
	 *            Array of row indices.
	 * @param j0
	 *            Initial column index
	 * @param j1
	 *            Final column index
	 * @param X
	 *            A(r(:),j0:j1)
	 * @exception ArrayIndexOutOfBoundsException
	 *                Submatrix indices
	 */

	public void setMatrix(int[] r, int j0, int j1, Matrix X) {
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = j0; j <= j1; j++) {
					_data[r[i] * _cols + j] = X.get(i, j - j0);
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
	}

	/**
	 * Set a submatrix.
	 * 
	 * @param i0
	 *            Initial row index
	 * @param i1
	 *            Final row index
	 * @param c
	 *            Array of column indices.
	 * @param X
	 *            A(i0:i1,c(:))
	 * @exception ArrayIndexOutOfBoundsException
	 *                Submatrix indices
	 */

	public void setMatrix(int i0, int i1, int[] c, Matrix X) {
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = 0; j < c.length; j++) {
					_data[i * _cols + c[j]] = X.get(i - i0, j);
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
	}

	/**
	 * Matrix rank
	 * 
	 * @return effective numerical rank, obtained from SVD.
	 */

	public int rank() {
		return new SingularValueDecomposition(this).rank();
	}

	/**
	 * Matrix condition (2 norm)
	 * 
	 * @return ratio of largest to smallest singular value.
	 */

	public double cond() {
		return new SingularValueDecomposition(this).cond();
	}

	/**
	 * Matrix trace.
	 * 
	 * @return sum of the diagonal elements.
	 */

	public double trace() {
		double t = 0;
		for (int i = 0; i < Math.min(_rows, _cols); ++i) {
			t += _data[i * _cols + i];
		}
		return t;
	}

	public Matrix cofactor() {
		int dim = this._rows;
		// Make sure that matrix is square

		double[][] cofactor = new double[dim][dim];

		int sign = -1;
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j < dim; ++j) {
				sign = -sign;
				cofactor[i][j] = sign * this.minor(i, j).det();
			}
		}
		return new Matrix(cofactor);
	}

	public boolean isSquared() {
		return _rows == _cols;
	}

	public Matrix minor(int row, int col) {
		final int dim = this._rows;
		double[][] minor = new double[dim - 1][dim - 1];

		for (int i = 0; i < dim; ++i) {
			int offset = i * this._cols;
			for (int j = 0; i != row && j < this._cols; ++j) {
				if (j != col) {
					minor[i < row ? i : i - 1][j < col ? j : j - 1] = this._data[offset + j];
				}
			}
		}
		return new Matrix(minor);
	}

	public Matrix adjoint() {
		return this.cofactor().transpose();
	}

	public LUDecomposition LU() {
		return new LUDecomposition(this);

	}

	public QRDecomposition QR() {
		return new QRDecomposition(this);
	}

	public CholeskyDecomposition Chol() {
		return new CholeskyDecomposition(this);
	}

	public SingularValueDecomposition SVD() {
		return new SingularValueDecomposition(this);
	}

	/***
	 * Get matrix diagonal
	 * 
	 * @return Array containing the diagonal of the matrix
	 */
	public double[] diag() {
		if (!this.isSquared()) {
			throw new RuntimeException("Matrix is not squared");
		}
		final int dim = _rows;
		double[] diag = new double[dim];
		for (int i = 0; i < dim; ++i) {
			diag[i] = _data[i * dim + i];
		}
		return diag;
	}

	public Matrix inv() {
		return this.solve(Matrices.Identity(_rows));
	}

	public Matrix transpose() {
		double[] result = new double[_rows * _cols];
		final int trows = _cols;
		final int tcols = _rows;

		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < _cols; ++j) {
				result[j * tcols + i] = _data[i * _cols + j];
			}
		}
		return new Matrix(result, trows, tcols);
	}

	/***
	 * One norm
	 * 
	 * @return maximum column sum
	 */
	public double norm1() {
		double norm = 0.0;
		for (int j = 0; j < _cols; ++j) {
			double sum = 0.0;
			for (int i = 0; i < _rows; ++i) {
				sum += _data[i * _cols + j];
			}
			norm = Math.max(norm, sum);
		}
		return norm;
	}

	/**
	 * Two norm
	 * 
	 * @return maximum singular value.
	 */

	public double norm2() {
		return new SingularValueDecomposition(this).norm2();
	}

	/***
	 * Infinity norm
	 * 
	 * @return maximum row sum
	 */
	public double normInf() {
		double norm = 0.0;
		for (int i = 0; i < _rows; ++i) {
			double sum = 0.0;
			for (int j = 0; j < _cols; ++j) {
				sum += _data[i * _cols + j];
			}
			norm = Math.max(norm, sum);
		}
		return norm;
	}

	/***
	 * Frobenius norm
	 * 
	 * @return square root of the sum of squares of all elements
	 */
	public double normFrob() {
		double norm = 0.0;
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < _cols; ++j) {
				norm = MathETK.hypot(norm, _data[i * _cols + j]);
			}
		}
		return norm;
	}

	/**
	 * Unary minus
	 * 
	 * @return -A
	 */
	public Matrix uminus() {
		final int m = _rows;
		final int n = _cols;
		final int length = m * n;
		double[] x = new double[length];
		for (int i = 0; i < length; ++i) {
			x[i] = -_data[i];
		}
		return new Matrix(x, m, n);
	}

	public Matrix add(Matrix m) {
		checkMatrixDimensions(m);
		double[] result = new double[this._rows * this._cols];
		for (int i = 0; i < this._rows * this._cols; ++i) {
			result[i] = this._data[i] + m._data[i];
		}
		return new Matrix(result, _rows, _cols);
	}

	public void addEquals(Matrix m) {
		checkMatrixDimensions(m);
		final int length = _rows * _cols;
		for (int i = 0; i < length; ++i) {
			_data[i] += m._data[i];
		}
	}

	public Matrix subtract(Matrix m) {
		checkMatrixDimensions(m);
		double[] result = new double[this._rows * this._cols];
		for (int i = 0; i < this._rows * this._cols; ++i) {
			result[i] = this._data[i] - m._data[i];
		}
		return new Matrix(result, _rows, _cols);
	}

	public void subtractEquals(Matrix m) {
		checkMatrixDimensions(m);
		final int length = _rows * _cols;
		for (int i = 0; i < length; ++i) {
			_data[i] -= m._data[i];
		}
	}

	/**
	 * Element-by-element multiplication, C = A.*B
	 * 
	 * @param B
	 *            another matrix
	 * @return A.*B
	 */

	public Matrix arrayMultiply(Matrix m) {
		checkMatrixDimensions(m);
		return new Matrix(NumArrays.multiplyElementWise(_data, m._data), _rows, _cols);
	}

	/**
	 * Element-by-element multiplication in place, A = A.*B
	 * 
	 * @param B
	 *            another matrix
	 * @return A.*B
	 */

	public void arrayMultiplyEquals(Matrix m) {
		checkMatrixDimensions(m);
		NumArrays.multiplyElementWiseInPlace(_data, m._data);
	}

	/**
	 * Element-by-element right division, C = A./B
	 * 
	 * @param B
	 *            another matrix
	 * @return A./B
	 */

	public Matrix arrayRightDivide(Matrix m) {
		checkMatrixDimensions(m);
		return new Matrix(NumArrays.divideElementWise(_data, m._data), _rows, _cols);
	}

	/**
	 * Element-by-element right division in place, A = A./B
	 * 
	 * @param B
	 *            another matrix
	 * @return A./B
	 */

	public void arrayRightDivideEquals(Matrix m) {
		checkMatrixDimensions(m);
		NumArrays.divideElementWiseInPlace(_data, m._data);
	}

	/**
	 * Element-by-element left division, C = A.\B
	 * 
	 * @param B
	 *            another matrix
	 * @return A.\B
	 */

	public Matrix arrayLeftDivide(Matrix m) {
		checkMatrixDimensions(m);
		return new Matrix(NumArrays.divideElementWise(m._data, _data), _rows, _cols);
	}

	/**
	 * Element-by-element left division in place, A = A.\B
	 * 
	 * @param B
	 *            another matrix
	 * @return A.\B
	 */

	public void arrayLeftDivideEquals(Matrix m) {
		checkMatrixDimensions(m);
		NumArrays.divideElementWiseInPlace(m._data, _data);
	}

	/**
	 * Multiply a matrix by a scalar, C = s*A
	 * 
	 * @param s
	 *            scalar
	 * @return s*A
	 */

	public Matrix multiply(double s) {
		return new Matrix(NumArrays.multiply(_data, s), _rows, _cols);
	}

	/**
	 * Multiply a matrix by a scalar in place, A = s*A
	 * 
	 * @param s
	 *            scalar
	 * @return replace A by s*A
	 */

	public void multiplyEquals(double s) {
		NumArrays.multiplyInPlace(_data, s);
	}

	public Matrix multiply(Matrix mat) {
		double[] result = new double[this._rows * mat._cols];
		for (int i = 0, j = 0; i < this._rows; i++) {
			for (int k = 0; k < mat._cols; k++, j++) {
				double value = this._data[i * this._cols] * mat._data[k];
				for (int l = 1, m = i * this._cols + 1, n = k + mat._cols; l < this._cols; l++, m++, n += mat._cols) {
					value += this._data[m] * mat._data[n];
				}
				result[j] = value;
			}
		}
		return new Matrix(result, _rows, mat._cols);
	}

	public Matrix solve(Matrix B) {

		if (_rows == _cols) { // Matrix is Squared
			return new LUDecomposition(this).solve(B);
		} else if (_rows > _cols) { // Matrix is thin (Overdetermined system)
			return new QRDecomposition(this).solve(B);
		} else { // Matrix is fat (Underdetermined system)
			QRDecomposition qr = this.transpose().QR();
			Matrix R1 = Matrices.fwdSubsSolve(qr.getRT(), B);
			R1.appendRows(_cols - R1._rows);
			return qr.QmultiplyX(R1);
		}
	}

	/**
	 * Solve X*A = B, which is also A'*X' = B'
	 * 
	 * @param B
	 *            right hand side
	 * @return solution if A is square, least squares solution otherwise.
	 */

	public Matrix solveTranspose(Matrix B) {
		return transpose().solve(B.transpose());
	}

	public void appendRows(int count) {
		_rows += count;
		final int newSize = _rows * _cols;
		_data = Arrays.copyOf(_data, newSize);

	}

	public Matrix pinv() {
		int rows = _rows;
		int cols = _cols;

		if (rows < cols) {
			Matrix result = this.transpose().pinv();
			if (result != null) {
				result = result.transpose();
			}
			return result;
		}

		SingularValueDecomposition svdX = this.SVD();
		if (svdX.rank() < 1) {
			return null;
		}

		double[] singularValues = svdX.getSingularValues();
		double tol = Math.max(rows, cols) * singularValues[0] * ETKConstants.DOUBLE_EPS;
		double[] singularValueReciprocals = new double[singularValues.length];
		for (int i = 0; i < singularValues.length; i++) {
			if (Math.abs(singularValues[i]) >= tol) {
				singularValueReciprocals[i] = 1.0 / singularValues[i];
			}
		}
		Matrix U = svdX.getU();
		Matrix V = svdX.getV();
		int min = Math.min(cols, U._cols);
		double[][] inverse = new double[cols][rows];
		for (int i = 0; i < cols; i++) {
			for (int j = 0; j < U._rows; j++) {
				for (int k = 0; k < min; k++) {
					inverse[i][j] += V.get(i, k) * singularValueReciprocals[k] * U.get(j, k);
				}
			}
		}
		return new Matrix(inverse);
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

	public EigenvalueDecomposition eig() {
		return this.eig(true);
	}

	public EigenvalueDecomposition eig(boolean balance) {
		Matrix A = balance ? new Matrix(this) : this;
		return new EigenvalueDecomposition(A, balance);
	}

	public int getRowCount() {
		return _rows;
	}

	public int getColumnCount() {
		return _cols;
	}

	public double[] getRow(int row) {
		double[] result = new double[_cols];
		int rowIndex = row * _cols;
		for (int j = 0; j < _cols; ++j) {
			result[j] = _data[rowIndex + j];
		}
		return result;
	}

	public double[] getArrayCopy() {
		return Arrays.copyOf(_data, _data.length);
	}

	public double[] getArray() {
		return _data;
	}

	public double[][] getAs2DArray() {
		double[][] data = new double[_rows][_cols];
		for (int i = 0; i < _rows; ++i) {
			data[i] = Arrays.copyOfRange(_data, i * _cols, i * _cols + _cols);
		}
		return data;
	}

	/**
	 * Make a one-dimensional column packed copy of the internal array.
	 * 
	 * @return Matrix elements packed in a one-dimensional array by columns.
	 */

	public double[] getColumnPackedCopy() {
		final int m = _rows;
		final int n = _cols;
		double[] vals = new double[m * n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				vals[i + j * m] = _data[i * n + j];
			}
		}
		return vals;
	}

	/**
	 * Make a one-dimensional row packed copy of the internal array.
	 * 
	 * @return Matrix elements packed in a one-dimensional array by rows.
	 */

	public double[] getRowPackedCopy() {
		return this.getArrayCopy();
	}

	/*
	 * ------------------------ Private Methods ------------------------
	 */

	/** Check if size(A) == size(B) **/

	private void checkMatrixDimensions(Matrix B) {
		if (B._rows != _rows || B._cols != _cols) {
			throw new IllegalArgumentException("Matrix dimensions must agree.");
		}
	}
}
