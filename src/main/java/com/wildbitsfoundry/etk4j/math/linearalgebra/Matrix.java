package com.wildbitsfoundry.etk4j.math.linearalgebra;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.curvefitting.CurveFitting;
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

    /***
     * Column packed
     * @param data
     * @param rows
     */
    public Matrix(double[] data, int rows) {
        _rows = rows;
        _cols = (_rows != 0 ? data.length / _rows : 0);
        if (_rows * _cols != data.length) {
            throw new IllegalArgumentException("Array length must be a multiple of rows");
        }

        int dim = _rows * _cols;
        _data = new double[dim];
        for (int i = 0; i < _rows; ++i) {
            for (int j = 0; j < _cols; ++j) {
                _data[i * _cols + j] = data[i + j * rows];
            }
        }
    }

    // row packed
    // _rows = rows
    // _cols = cols
    // _data = data;

    public Matrix(double[][] data) {
        _rows = data.length;
        _cols = data[0].length;
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
     * Deep copy
     * @return
     */
    public Matrix copy() {
        double[] data = Arrays.copyOf(this._data, this._data.length);
        return new Matrix(data, this._rows, this._cols);
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
        // check bounds
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
     * @return A(rows ( :), col0 : col1)
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
     * @return A(row0 : row1, cols ( :))
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
     * @return A(rows ( :), cols(:))
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
     * @param i0 Initial row index
     * @param i1 Final row index
     * @param j0 Initial column index
     * @param j1 Final column index
     * @param X  A(i0:i1,j0:j1)
     * @throws ArrayIndexOutOfBoundsException Submatrix indices
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
     * @param r Array of row indices.
     * @param c Array of column indices.
     * @param X A(r(:),c(:))
     * @throws ArrayIndexOutOfBoundsException Submatrix indices
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
     * @param r  Array of row indices.
     * @param j0 Initial column index
     * @param j1 Final column index
     * @param X  A(r(:),j0:j1)
     * @throws ArrayIndexOutOfBoundsException Submatrix indices
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
     * @param i0 Initial row index
     * @param i1 Final row index
     * @param c  Array of column indices.
     * @param X  A(i0:i1,c(:))
     * @throws ArrayIndexOutOfBoundsException Submatrix indices
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
                sum += Math.abs(_data[i * _cols + j]);
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
                sum += Math.abs(_data[i * _cols + j]);
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

    // TODO documents all this multiply methods. Change m for B
    /**
     * Element-by-element multiplication, C = A.*B
     *
     * @param m another matrix
     * @return A.*B
     */

    public Matrix arrayMultiply(Matrix m) {
        checkMatrixDimensions(m);
        return new Matrix(NumArrays.multiplyElementWise(_data, m._data), _rows, _cols);
    }

    /**
     * Element-by-element multiplication in place, A = A.*B
     *
     * @param m another matrix
     * @return A.*B
     */

    public void arrayMultiplyEquals(Matrix m) {
        checkMatrixDimensions(m);
        NumArrays.multiplyElementWiseInPlace(_data, m._data);
    }

    /**
     * Element-by-element right division, C = A./B
     *
     * @param m another matrix
     * @return A./B
     */

    public Matrix arrayRightDivide(Matrix m) {
        checkMatrixDimensions(m);
        return new Matrix(NumArrays.divideElementWise(_data, m._data), _rows, _cols);
    }

    /**
     * Element-by-element right division in place, A = A./B
     *
     * @param m another matrix
     * @return A./B
     */

    public void arrayRightDivideEquals(Matrix m) {
        checkMatrixDimensions(m);
        NumArrays.divideElementWiseInPlace(_data, m._data);
    }

    /**
     * Element-by-element left division, C = A.\B
     *
     * @param m another matrix
     * @return A.\B
     */

    public Matrix arrayLeftDivide(Matrix m) {
        checkMatrixDimensions(m);
        if (_data.length != m._data.length) {
            throw new IllegalArgumentException("a and b dimensions must match");
        }

        final int length = _data.length;
        double[] data = Arrays.copyOf(_data, length);
        for (int i = 0; i < length; ++i) {
            data[i] = m._data[i] / data[i];
        }
        return new Matrix(data, _rows, _cols);
    }

    /**
     * Element-by-element left division in place, A = A.\B
     *
     * @param m another matrix
     * @return A.\B
     */

    public void arrayLeftDivideEquals(Matrix m) {
        checkMatrixDimensions(m);
        if (_data.length != m._data.length) {
            throw new IllegalArgumentException("a and b dimensions must match");
        }

        final int length = _data.length;
        for (int i = 0; i < length; ++i) {
            _data[i] = m._data[i] / _data[i];
        }
    }

    /**
     * Multiply a matrix by a scalar, C = s*A
     *
     * @param s scalar
     * @return s*A
     */

    public Matrix multiply(double s) {
        return new Matrix(NumArrays.multiply(_data, s), _rows, _cols);
    }

    /**
     * Multiply a matrix by a scalar in place, A = s*A
     *
     * @param s scalar
     * @return replace A by s*A
     */

    public void multiplyEquals(double s) {
        NumArrays.multiplyInPlace(_data, s);
    }

    public Matrix multiply(Matrix matrix) {
        double[] data = null;
        Matrix c = new Matrix(0, 0);
        multiplyOp(this, matrix, c);
        return c;
    }

    private static void multiplyOp(Matrix a, Matrix b, Matrix c) {
        if (b._rows != a._cols) {
            throw new IllegalArgumentException("Matrix inner dimensions must agree.");
        }
        double[] result = new double[a._rows * b._cols];
        double[] bColJ = new double[a._cols];
        for (int j = 0; j < b._cols; j++) {
            for (int k = 0; k < a._cols; k++) {
                bColJ[k] = b.get(k, j);
            }
            for (int i = 0; i < a._rows; i++) {
                double[] aRowI = new double[a._cols];
                System.arraycopy(a._data, i * a._cols, aRowI, 0, a._cols);

                double s = 0;
                for (int k = 0; k < a._cols; k++) {
                    s += aRowI[k] * bColJ[k];
                }
                result[i * b._cols + j] = s;
            }
        }
        c._data = result;
        c._rows = a._rows;
        c._cols = b._cols;
    }

    public void multiplyEquals(Matrix matrix) {
        multiplyOp(this, matrix, this);
    }

    public Matrix solve(Matrix B) {

        if (_rows == _cols) { // Matrix is Squared
            return new LUDecomposition(this).solve(B);
        } else if (_rows > _cols) { // Matrix is thin (Overdetermined system)
            return new QRDecomposition(this).solve(B);
        } else { // Matrix is fat (Underdetermined system)
            QRDecomposition qr = this.transpose().QR();
            Matrix R1 = fwdSubsSolve(qr.getRT(), B);
            R1.appendRows(_cols - R1._rows);
            return qr.QmultiplyX(R1);
        }
    }

    private static Matrix fwdSubsSolve(Matrix L, Matrix B) {
        final int nx = B.getColumnCount();
        final int m = L.getRowCount();
        final int n = L.getColumnCount();
        double[] t = L.getArray();
        double[] X = B.getArrayCopy();

        for (int j = 0; j < nx; ++j) {
            for (int i = 0; i < m; ++i) {
                for (int k = 0; k < i; ++k) {
                    X[i * nx + j] -= X[k * nx + j] * t[i * n + k];
                }
                X[i * nx + j] /= t[i * n + i];
            }
        }
        return new Matrix(X, m, nx);
    }

    /**
     * Solve X*A = B, which is also A'*X' = B'
     *
     * @param B right hand side
     * @return solution if A is square, least squares solution otherwise.
     */

    public Matrix transposeSolve(Matrix B) {
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
        double tol = Math.max(rows, cols) * singularValues[0] * ConstantsETK.DOUBLE_EPS;
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
            sb.append(String.format("%.4g", _data[i])).append(" ");
        }
        sb.setLength(sb.length() - 1);
        return sb.toString();
    }

    public EigenvalueDecomposition eig() {
        return new EigenvalueDecomposition(this);
    }

    public EigenvalueDecomposition eig(boolean balance) {
        return new EigenvalueDecomposition(this, balance);
    }

    /***
     * Balances the matrix using the algorithm by parlett and reinsch with norm - 1.
     * http://www.netlib.org/eispack/balanc.f
     * https://arxiv.org/pdf/1401.5766.pdf algorithm #2
     *
     * @return
     *
     * @see "https://stackoverflow.com/a/43169781/6383857"
     */
    public Matrix balance() {
        if (!this.isSquared()) {
            throw new IllegalArgumentException("Matrix must be squared.");
        }
        int rows = _rows;
        double[] data = this.getArrayCopy();
        double radix = 2.0;                // radix base
        double radix2 = radix * radix;    // radix base squared
        boolean done = false;
        while (!done) {
            done = true;
            for (int i = 0; i < rows; i++) {
                double r = 0.0, c = 0.0;
                for (int j = 0; j < rows; j++) {
                    if (j != i) {
                        // Compute row[i] and col[i] norm - 1
                        c += Math.abs(data[j * rows + i]);
                        r += Math.abs(data[i * rows + j]);
                    }
                }
                if (c != 0 && r != 0) {
                    double s = c + r;
                    double f = 1.0;
                    double g = r / radix;
                    while (c < g) {
                        f *= radix;
                        c *= radix2;
                    }
                    g = r * radix;
                    while (c > g) {
                        f /= radix;
                        c /= radix2;
                    }
                    if ((c + r) / f < 0.95 * s) {
                        done = false;
                        g = 1.0 / f;
                        //scaling[i] *= f;
                        for (int j = 0; j < rows; j++) {
                            data[i * rows + j] *= g;
                        }
                        for (int j = 0; j < rows; j++) {
                            data[j * rows + i] *= f;
                        }
                    }
                }
            }
        }
//        double beta = 2;
//        int rows = this.getRowCount();
//        double[] data = this.getArrayCopy();
//        boolean converged = false;
//        do {
//            converged = true;
//            for(int i = 0; i < rows; ++i) {
//                double c = 0.0;
//                double r = 0.0;
//                // norm - 2
//                for (int j = 0; j < rows; ++j) {
//                    c += data[j * rows + i] * data[j * rows + i];
//                    r += data[i * rows + j] * data[i * rows + j];
//                }
//                c = Math.sqrt(c);
//                r = Math.sqrt(r);
//                if(c != 0.0 && r != 0.0) {
//                    double f = 1.0;
//                    double s = c * c + r * r;
//                    double alpha = 1.0 / beta;
//                    while (c < r / beta) {
//                        c *= beta;
//                        r *= alpha;
//                        f *= beta;
//                    }
//                    while (c >= r * beta) {
//                        c *= alpha;
//                        r *= beta;
//                        f *= alpha;
//                    }
//                    if (c * c + r * r < 0.95 * s) {
//                        converged = false;
//                        double g = 1.0 / f;
//                        for (int j = 0; j < rows; ++j) {
//                            data[j * rows + i] *= f;
//                            data[i * rows + j] *= g;
//                        }
//                    }
//                }
//            }
//        } while(!converged);
        return new Matrix(data, rows, rows);
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

    public double[] getCol(int col) {
        double[] result = new double[_rows];
        for (int i = 0; i < _rows; ++i) {
            result[i] = _data[i * _rows + col];
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


    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + _cols;
        result = prime * result + Arrays.hashCode(_data);
        result = prime * result + _rows;
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (!(obj instanceof Matrix)) {
            return false;
        }
        Matrix other = (Matrix) obj;
        if (_cols != other._cols) {
            return false;
        }
        if (!Arrays.equals(_data, other._data)) {
            return false;
        }
        if (_rows != other._rows) {
            return false;
        }
        return true;
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

    /**
     * Check if size(A) == size(B)
     **/

    private void checkMatrixDimensions(Matrix B) {
        if (B._rows != _rows || B._cols != _cols) {
            throw new IllegalArgumentException("Matrix dimensions must agree.");
        }
    }

    public void setRow(int i, double[] row) {
        if (i < 0 || i >= _rows) {
            throw new IndexOutOfBoundsException("f");
        }
        if (row.length != _cols) {
            throw new IllegalArgumentException();
        }
        System.arraycopy(row, 0, _data, i * _cols, _cols);
    }

    public Matrix pow(int n) {
        if (!this.isSquared()) {
            throw new IllegalArgumentException("Matrix must be squared");
        }
        n = Math.abs(n);
        if (n == 0) {
            return Matrices.Identity(_rows);
        }
        Matrix a = n < 0 ? inv() : new Matrix(this);

        if (n == 1) {
            return a;
        } else if (n == 2) {
            a.multiplyEquals(a);
            return a;
        } else if (n == 3) {
            a.multiplyEquals(a.multiply(a));
            return a;
        }
        Matrix z = null, result = null;
        while(n > 0) {
            if(z == null) {
                z = a;
            } else {
                z = z.multiply(z);
            }
            int bit = n % 2;
            n = n / 2;
            if(bit > 0) {
                if(result == null) {
                    result = z;
                } else {
                    result.multiplyEquals(z);
                }
            }
        }
        return result;
    }

    public Matrix expm() {
        FRexpResult result = frexp(this.normInf());
        double s = Math.max(0.0, result.exponent + 1);
        Matrix A = this.multiply(1 / Math.pow(2.0, s));

        Matrix X = new Matrix(A);
        double c = 0.5;
        Matrix E = Matrices.Identity(A._rows, A._cols).add(A.multiply(c));
        Matrix D = Matrices.Identity(A._rows, A._cols).subtract(A.multiply(c));
        double q = 6.0;
        boolean p = true;
        for(int k = 2; k <= q; ++k) {
            c = c * (q - k + 1) / (k * (2 * q - k + 1));
            X.multiplyEquals(A);
            Matrix cX = X.multiply(c);
            E.addEquals(cX);
            if(p) {
                D.addEquals(cX);
            } else {
                D.subtractEquals(cX);
            }
            p = !p;
        }
        E = D.solve(E);

        for(int k = 0; k < s; ++k) {
            E.multiplyEquals(E);
        }
        return E;

    }

//    public Matrix pow(double n) {
//        if (!this.isSquared()) {
//            throw new IllegalArgumentException("Matrix must be squared");
//        }
//        n = Math.abs(n);
//        if (n == 0.0 || n == 1.0 || n == 2.0 || n == 3.0) {
//            return pow((int) n);
//        }
//        Matrix a = n < 0 ? inv() : new Matrix(this);
//        EigenvalueDecomposition eig = a.eig();
//        Matrix D = eig.getD();
//        for(int i = 0, j = 0; i < D._rows; ++i, ++j) {
//            D.set(i, j, Math.pow(D.get(i, j), n));
//        }
//        Matrix result = eig.getV().multiply(D);
//        result.multiplyEquals(eig.getV().inv());
//        return result;
//    }

    // TODO add more matrix multiplication tests?
    public static void main(String[] args) {
        Matrix a = Matrices.Magic(3);
        double[] rowPacked = a.getRowPackedCopy();
        double[] colPacked = a.getColumnPackedCopy();

        System.out.println(a);
        System.out.println(Arrays.toString(rowPacked));
        System.out.println(Arrays.toString(colPacked));

        System.out.println(Arrays.toString(a.getRow(1)));
        System.out.println(Arrays.toString(a.getCol(2)));

        Matrix b = Matrices.Magic(3);
        Matrix c = new Matrix(b);

        b = b.balance();

        System.out.println();
        System.out.println(b);
        System.out.println();
        System.out.println(c);

        double[] data = {1, 0.01, 0.0001, 100, 1, 0.01, 10000, 100, 1};
        Matrix d = new Matrix(data, 3);
        System.out.println();
        System.out.println(d);

        d = d.balance();
        System.out.println();
        System.out.println(d);

        Matrix mul = Matrices.Magic(3);
        System.out.println();
        System.out.println(mul);
        System.out.println();
        System.out.println(mul.pow(0));
        System.out.println();
        System.out.println(mul.pow(1));
        System.out.println();
        System.out.println(mul.pow(2));
        System.out.println();
        System.out.println(mul.pow(3));
        System.out.println();
        System.out.println(mul.pow(4));
        System.out.println();
        System.out.println(mul.pow(5));

        double[][] ddata = {{1,2,3}, {4,5,6},{7,8,9}};

        Matrix E = new Matrix(ddata);
        System.out.println();
        System.out.println(E.expm());
    }

    public static class FRexpResult
    {
        public int exponent = 0;
        public double mantissa = 0.;
    }

    public static FRexpResult frexp(double value)
    {
        final FRexpResult result = new FRexpResult();
        long bits = Double.doubleToLongBits(value);
        double realMant = 1.;

        // Test for NaN, infinity, and zero.
        if (Double.isNaN(value) ||
                value + value == value ||
                Double.isInfinite(value))
        {
            result.exponent = 0;
            result.mantissa = value;
        }
        else
        {

            boolean neg = (bits < 0);
            int exponent = (int)((bits >> 52) & 0x7ffL);
            long mantissa = bits & 0xfffffffffffffL;

            if(exponent == 0)
            {
                exponent++;
            }
            else
            {
                mantissa = mantissa | (1L<<52);
            }

            // bias the exponent - actually biased by 1023.
            // we are treating the mantissa as m.0 instead of 0.m
            //  so subtract another 52.
            exponent -= 1075;
            realMant = mantissa;

            // normalize
            while(realMant >= 1.0)
            {
                mantissa >>= 1;
                realMant /= 2.;
                exponent++;
            }

            if(neg)
            {
                realMant = realMant * -1;
            }

            result.exponent = exponent;
            result.mantissa = realMant;
        }
        return result;
    }
}
