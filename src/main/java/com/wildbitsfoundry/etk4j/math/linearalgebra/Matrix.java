package com.wildbitsfoundry.etk4j.math.linearalgebra;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.NumArrays;

import static com.wildbitsfoundry.etk4j.math.MathETK.frexp;

public class Matrix {
    private double[] data;
    private int rows;
    private int cols;

    public Matrix(int rows, int cols) {
        this.rows = rows;
        this.cols = cols;

        this.data = new double[rows * cols];
    }

    /***
     * Column packed
     * @param data
     * @param rows
     */
    public Matrix(double[] data, int rows) {
        this.rows = rows;
        cols = (this.rows != 0 ? data.length / this.rows : 0);
        if (this.rows * cols != data.length) {
            throw new IllegalArgumentException("Array length must be a multiple of rows");
        }

        int dim = this.rows * cols;
        this.data = new double[dim];
        for (int i = 0; i < this.rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                this.data[i * cols + j] = data[i + j * rows];
            }
        }
    }

    public Matrix(double[][] data) {
        rows = data.length;
        cols = data[0].length;
        this.data = NumArrays.flatten(data);
    }

    // row packed
    // _rows = rows
    // _cols = cols
    // _data = data;
    public Matrix(double[] data, int rows, int cols) {
        this.rows = rows;
        this.cols = cols;
        this.data = data;
    }

    public Matrix(Matrix matrix) {
        rows = matrix.rows;
        cols = matrix.cols;
        data = new double[rows * cols];
        System.arraycopy(matrix.data, 0, this.data, 0, this.rows * this.cols);
    }

    public Matrix(int rows, int cols, double val) {
        this.rows = rows;
        this.cols = cols;
        data = new double[this.rows * this.cols];
        Arrays.fill(data, val);
    }

    /***
     * Deep copy
     * @return
     */
    public Matrix copy() {
        double[] data = Arrays.copyOf(this.data, this.data.length);
        return new Matrix(data, this.rows, this.cols);
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
        // TODO check bounds and check that the 1's are bigger than the 0's
        // check bounds
        int rowDim = row1 - row0 + 1;
        int colDim = col1 - col0 + 1;
        double[] data = new double[rowDim * colDim];
        for (int i = row0; i <= row1; ++i) {
            for (int j = col0; j <= col1; ++j) {
                data[(i - row0) * colDim + (j - col0)] = this.data[i * cols + j];
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
                data[i * colDim + j] = this.data[rows[i] * cols + (j + col0)];
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
                data[i * colDim + j] = this.data[(i + row0) * this.cols + cols[j]];
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
                data[i * colDim + j] = this.data[rows[i] * this.cols + cols[j]];
            }
        }
        return new Matrix(data, rowDim, colDim);
    }

    public double get(int i, int j) {
        if(i < 0) {
            throw new ArrayIndexOutOfBoundsException("Index i cannot be less thant zero.");
        }
        if (i >= rows) {
            throw new ArrayIndexOutOfBoundsException(String.format("Index i: %d >= than number of rows: %d.", i, rows));
        }
        if(j < 0) {
            throw new ArrayIndexOutOfBoundsException("Index j cannot be less thant zero.");
        }
        if (j >= cols) {
            throw new ArrayIndexOutOfBoundsException(String.format("Index j: %d >= than number of columns: %d.", j, cols));
        }
        return data[i * cols + j];
    }

    public void set(int i, int j, double val) {
        if(i < 0) {
            throw new ArrayIndexOutOfBoundsException("Index i cannot be less thant zero.");
        }
        if (i >= rows) {
            throw new ArrayIndexOutOfBoundsException(String.format("Index i: %d >= than number of rows: %d.", i, rows));
        }
        if(j < 0) {
            throw new ArrayIndexOutOfBoundsException("Index j cannot be less thant zero.");
        }
        if (j >= cols) {
            throw new ArrayIndexOutOfBoundsException(String.format("Index j: %d >= than number of columns: %d.", j, cols));
        }
        data[i * cols + j] = val;
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
                    data[i * cols + j] = X.get(i - i0, j - j0);
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
                    data[r[i] * cols + c[j]] = X.get(i, j);
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
                    data[r[i] * cols + j] = X.get(i, j - j0);
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
                    data[i * cols + c[j]] = X.get(i - i0, j);
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
        for (int i = 0; i < Math.min(rows, cols); ++i) {
            t += data[i * cols + i];
        }
        return t;
    }

    public Matrix cofactor() {
        int dim = this.rows;
        if(!this.isSquared()) {
            throw new NonSquareMatrixException("Matrix must be a square Matrix.");
        }
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
        return rows == cols;
    }

    public Matrix minor(int row, int col) {
        final int dim = this.rows;
        double[][] minor = new double[dim - 1][dim - 1];

        for (int i = 0; i < dim; ++i) {
            int offset = i * this.cols;
            for (int j = 0; i != row && j < this.cols; ++j) {
                if (j != col) {
                    minor[i < row ? i : i - 1][j < col ? j : j - 1] = this.data[offset + j];
                }
            }
        }
        return new Matrix(minor);
    }

    public Matrix adjoint() {
        return this.cofactor().transpose();
    }

    // region decompositions
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
    // endregion
    /***
     * Get matrix diagonal
     *
     * @return Array containing the diagonal of the matrix
     */
    public double[] diag() {
        if (!this.isSquared()) {
            throw new RuntimeException("Matrix is not squared");
        }
        final int dim = rows;
        double[] diag = new double[dim];
        for (int i = 0; i < dim; ++i) {
            diag[i] = data[i * dim + i];
        }
        return diag;
    }

    public Matrix inv() {
        return this.solve(Matrices.identity(rows));
    }

    public boolean isEmpty() {
        if ((rows == 0 && cols == 0) || data == null || data.length == 0) {
            return true;
        }
        return false;
    }

    public Matrix transpose() {
        if (this.isEmpty()) {
            return Matrices.empty();
        }
        double[] result = new double[rows * cols];
        final int trows = cols;
        final int tcols = rows;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result[j * tcols + i] = data[i * cols + j];
            }
        }
        return new Matrix(result, trows, tcols);
    }

    // region norms
    /***
     * One norm
     *
     * @return maximum column sum
     */
    public double norm1() {
        double norm = 0.0;
        for (int j = 0; j < cols; ++j) {
            double sum = 0.0;
            for (int i = 0; i < rows; ++i) {
                sum += Math.abs(data[i * cols + j]);
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
        for (int i = 0; i < rows; ++i) {
            double sum = 0.0;
            for (int j = 0; j < cols; ++j) {
                sum += Math.abs(data[i * cols + j]);
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
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                norm = MathETK.hypot(norm, data[i * cols + j]);
            }
        }
        return norm;
    }
    // endregion

    /**
     * Unary minus
     *
     * @return -A
     */
    public Matrix uminus() {
        final int m = rows;
        final int n = cols;
        final int length = m * n;
        double[] x = new double[length];
        for (int i = 0; i < length; ++i) {
            x[i] = -data[i];
        }
        return new Matrix(x, m, n);
    }

    // region arithmetic operations
    public Matrix add(Matrix m) {
        checkMatrixDimensions(m);
        double[] result = new double[this.rows * this.cols];
        for (int i = 0; i < this.rows * this.cols; ++i) {
            result[i] = this.data[i] + m.data[i];
        }
        return new Matrix(result, rows, cols);
    }

    public void addEquals(Matrix m) {
        checkMatrixDimensions(m);
        final int length = rows * cols;
        for (int i = 0; i < length; ++i) {
            data[i] += m.data[i];
        }
    }

    public Matrix subtract(Matrix m) {
        checkMatrixDimensions(m);
        double[] result = new double[this.rows * this.cols];
        for (int i = 0; i < this.rows * this.cols; ++i) {
            result[i] = this.data[i] - m.data[i];
        }
        return new Matrix(result, rows, cols);
    }

    public void subtractEquals(Matrix m) {
        checkMatrixDimensions(m);
        final int length = rows * cols;
        for (int i = 0; i < length; ++i) {
            data[i] -= m.data[i];
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
        return new Matrix(NumArrays.multiplyElementWise(data, m.data), rows, cols);
    }

    /**
     * Element-by-element multiplication in place, A = A.*B
     *
     * @param m another matrix
     * @return A.*B
     */

    public void arrayMultiplyEquals(Matrix m) {
        checkMatrixDimensions(m);
        NumArrays.multiplyElementWiseInPlace(data, m.data);
    }

    /**
     * Element-by-element right division, C = A./B
     *
     * @param m another matrix
     * @return A./B
     */

    public Matrix arrayRightDivide(Matrix m) {
        checkMatrixDimensions(m);
        return new Matrix(NumArrays.divideElementWise(data, m.data), rows, cols);
    }

    /**
     * Element-by-element right division in place, A = A./B
     *
     * @param m another matrix
     * @return A./B
     */

    public void arrayRightDivideEquals(Matrix m) {
        checkMatrixDimensions(m);
        NumArrays.divideElementWiseInPlace(data, m.data);
    }

    /**
     * Element-by-element left division, C = A.\B
     *
     * @param m another matrix
     * @return A.\B
     */

    public Matrix arrayLeftDivide(Matrix m) {
        checkMatrixDimensions(m);
        if (data.length != m.data.length) {
            throw new IllegalArgumentException("a and b dimensions must match");
        }

        final int length = data.length;
        double[] data = Arrays.copyOf(this.data, length);
        for (int i = 0; i < length; ++i) {
            data[i] = m.data[i] / data[i];
        }
        return new Matrix(data, rows, cols);
    }

    /**
     * Element-by-element left division in place, A = A.\B
     *
     * @param m another matrix
     * @return A.\B
     */

    public void arrayLeftDivideEquals(Matrix m) {
        checkMatrixDimensions(m);
        if (data.length != m.data.length) {
            throw new IllegalArgumentException("a and b dimensions must match");
        }

        final int length = data.length;
        for (int i = 0; i < length; ++i) {
            data[i] = m.data[i] / data[i];
        }
    }

    /**
     * Multiply a matrix by a scalar, C = s*A
     *
     * @param s scalar
     * @return s*A
     */

    public Matrix multiply(double s) {
        return new Matrix(NumArrays.multiplyElementWise(data, s), rows, cols);
    }

    /**
     * Multiply a matrix by a complex scalar, C = s*A
     *
     * @param s scalar
     * @return s*A
     */

    public ComplexMatrix multiply(Complex s) {
        return new ComplexMatrix(ComplexArrays.multiplyElementWise(data, s), rows, cols);
    }

    /**
     * Multiply a matrix by a scalar in place, A = s*A
     *
     * @param s scalar
     * @return replace A by s*A
     */

    public void multiplyEquals(double s) {
        NumArrays.multiplyElementWiseInPlace(data, s);
    }

    public Matrix multiply(Matrix matrix) {
        Matrix c = new Matrix(0, 0);
        multiplyOp(this, matrix, c);
        return c;
    }

    public void multiplyEquals(Matrix matrix) {
        multiplyOp(this, matrix, this);
    }

    private static void multiplyOp(Matrix a, Matrix b, Matrix c) {
        if (b.rows != a.cols) {
            throw new IllegalArgumentException("Matrix inner dimensions must agree. Check that the number of" +
                    "columns of the first matrix equal the number of rows of the second matrix.");
        }
        double[] result = new double[a.rows * b.cols];
        double[] bColJ = new double[a.cols];
        for (int j = 0; j < b.cols; j++) {
            for (int k = 0; k < a.cols; k++) {
                bColJ[k] = b.data[k * b.cols + j];
            }
            for (int i = 0; i < a.rows; i++) {
                double s = 0;
                for (int k = 0; k < a.cols; k++) {
                    s += a.data[k + i * a.cols] * bColJ[k];
                }
                result[i * b.cols + j] = s;
            }
        }
        c.data = result;
        c.rows = a.rows;
        c.cols = b.cols;
    }

    public ComplexMatrix multiply(ComplexMatrix matrix) {
        int bRows = matrix.getRowCount();
        int bCols = matrix.getColumnCount();
        if (bRows != cols) {
            throw new IllegalArgumentException("Matrix inner dimensions must agree. Check that the number of" +
                    "columns of the first matrix equal the number of rows of the second matrix.");
        }
        Complex[] result = new Complex[rows * bCols];
        Complex[] bColJ = new Complex[cols];
        Complex[] bData = matrix.getArray();
        for (int j = 0; j < bCols; j++) {
            for (int k = 0; k < cols; k++) {
                bColJ[k] = bData[k * bCols + j];
            }
            for (int i = 0; i < rows; i++) {
                Complex s = new Complex();
                for (int k = 0; k < cols; k++) {
                    s.addEquals(bColJ[k].multiply(data[k + i * cols]));
                }
                result[i * bCols + j] = s;
            }
        }
        return new ComplexMatrix(result, rows, bCols);
    }
    // endregion
    public Matrix solve(Matrix B) {

        if (rows == cols) { // Matrix is Squared
            return new LUDecomposition(this).solve(B);
        } else if (rows > cols) { // Matrix is thin (Overdetermined system)
            return new QRDecomposition(this).solve(B);
        } else { // Matrix is fat (Under-determined system)
            QRDecomposition qr = this.transpose().QR();
            Matrix R1 = fwdSubsSolve(qr.getRT(), B);
            R1.appendRows(cols - R1.rows);
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
        rows += count;
        final int newSize = rows * cols;
        data = Arrays.copyOf(data, newSize);

    }

    public Matrix pinv() {
        int rows = this.rows;
        int cols = this.cols;

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
        int min = Math.min(cols, U.cols);
        double[][] inverse = new double[cols][rows];
        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < U.rows; j++) {
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
        for (int i = 0; i < rows * cols; ++i) {
            if (i > 0 && i % cols == 0) {
                sb.append(System.lineSeparator());
            }
            sb.append(String.format("%.4g", data[i])).append(" ");
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
            throw new NonSquareMatrixException("Matrix must be a square Matrix.");
        }
        int rows = this.rows;
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
        return rows;
    }

    public int getColumnCount() {
        return cols;
    }

    public double[] getRow(int row) {
        double[] result = new double[cols];
        int rowIndex = row * cols;
        for (int j = 0; j < cols; ++j) {
            result[j] = data[rowIndex + j];
        }
        return result;
    }

    public double[] getCol(int col) {
        double[] result = new double[rows];
        for (int i = 0; i < rows; ++i) {
            result[i] = data[i * rows + col];
        }
        return result;
    }

    public double[] getArrayCopy() {
        return Arrays.copyOf(data, data.length);
    }

    public double[] getArray() {
        return data;
    }

    public double[][] getAs2DArray() {
        double[][] data = new double[rows][cols];
        for (int i = 0; i < rows; ++i) {
            data[i] = Arrays.copyOfRange(this.data, i * cols, i * cols + cols);
        }
        return data;
    }

    /**
     * Make a one-dimensional column packed copy of the internal array.
     *
     * @return Matrix elements packed in a one-dimensional array by columns.
     */

    public double[] getColumnPackedCopy() {
        final int m = rows;
        final int n = cols;
        double[] vals = new double[m * n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                vals[i + j * m] = data[i * n + j];
            }
        }
        return vals;
    }


    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + cols;
        result = prime * result + Arrays.hashCode(data);
        result = prime * result + rows;
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
        if (cols != other.cols) {
            return false;
        }
        if (!Arrays.equals(data, other.data)) {
            return false;
        }
        if (rows != other.rows) {
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
    // TODO should this be static?
    void checkMatrixDimensions(Matrix B) {
        if (B.rows != rows || B.cols != cols) {
            throw new IllegalArgumentException("Matrix dimensions must agree.");
        }
    }

    // TODO
    public void setRow(int i, double[] row) {
        if (i < 0 || i >= rows) {
            throw new IndexOutOfBoundsException("f");
        }
        if (row.length != cols) {
            throw new IllegalArgumentException();
        }
        System.arraycopy(row, 0, data, i * cols, cols);
    }

    /**
     * Power of Matrix.
     * @param n The power to raise the matrix to.
     * @return A^n
     */
    public Matrix pow(int n) {
        if (!this.isSquared()) {
            throw new IllegalArgumentException("Matrix must be a square matrix.");
        }
        n = Math.abs(n);
        if (n == 0) {
            return Matrices.identity(rows);
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
        while (n > 0) {
            if (z == null) {
                z = a;
            } else {
                z = z.multiply(z);
            }
            int bit = n % 2;
            n = n / 2;
            if (bit > 0) {
                if (result == null) {
                    result = z;
                } else {
                    result.multiplyEquals(z);
                }
            }
        }
        return result;
    }

    /**
     * Exponential of Matrix.
     * @return e^A
     */
    public Matrix expm() {
        MathETK.FRexpResult result = frexp(this.normInf());
        double s = Math.max(0.0, result.exponent + 1);
        Matrix A = this.multiply(1 / Math.pow(2.0, s));

        Matrix X = new Matrix(A);
        double c = 0.5;
        Matrix E = Matrices.identity(A.rows, A.cols).add(A.multiply(c));
        Matrix D = Matrices.identity(A.rows, A.cols).subtract(A.multiply(c));
        double q = 6.0;
        boolean p = true;
        for (int k = 2; k <= q; ++k) {
            c = c * (q - k + 1) / (k * (2 * q - k + 1));
            X.multiplyEquals(A);
            Matrix cX = X.multiply(c);
            E.addEquals(cX);
            if (p) {
                D.addEquals(cX);
            } else {
                D.subtractEquals(cX);
            }
            p = !p;
        }
        E = D.solve(E);

        for (int k = 0; k < s; ++k) {
            E.multiplyEquals(E);
        }
        return E;

    }

    // Power using Eigen values
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
//        Matrix a = Matrices.Magic(3);
//        double[] rowPacked = a.getRowPackedCopy();
//        double[] colPacked = a.getColumnPackedCopy();
//
//        System.out.println(a);
//        System.out.println(Arrays.toString(rowPacked));
//        System.out.println(Arrays.toString(colPacked));
//
//        System.out.println(Arrays.toString(a.getRow(1)));
//        System.out.println(Arrays.toString(a.getCol(2)));
//
//        Matrix b = Matrices.Magic(3);
//        Matrix c = new Matrix(b);
//
//        b = b.balance();
//
//        System.out.println();
//        System.out.println(b);
//        System.out.println();
//        System.out.println(c);
//
//        double[] data = {1, 0.01, 0.0001, 100, 1, 0.01, 10000, 100, 1};
//        Matrix d = new Matrix(data, 3);
//        System.out.println();
//        System.out.println(d);
//
//        d = d.balance();
//        System.out.println();
//        System.out.println(d);
//
//        Matrix mul = Matrices.Magic(3);
//        System.out.println();
//        System.out.println(mul);
//        System.out.println();
//        System.out.println(mul.pow(0));
//        System.out.println();
//        System.out.println(mul.pow(1));
//        System.out.println();
//        System.out.println(mul.pow(2));
//        System.out.println();
//        System.out.println(mul.pow(3));
//        System.out.println();
//        System.out.println(mul.pow(4));
//        System.out.println();
//        System.out.println(mul.pow(5));
//
//        double[][] ddata = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//
//        Matrix E = new Matrix(ddata);
//        System.out.println();
//        System.out.println(E.expm());
//
//        double[][] Aa = {{-2, -1}, {1, 0}};
//        EigenvalueDecomposition eigenvalueDecomposition = new EigenvalueDecomposition(new Matrix(Aa));
//        Complex[] roots = ComplexArrays.zip(eigenvalueDecomposition.getRealEigenvalues(), eigenvalueDecomposition.getImagEigenvalues());
//        double[] coefficients = new Polynomial(roots).getCoefficients();
//        System.out.println(Arrays.toString(coefficients));
//
//        Aa = new double[][]{{-1}};
//        eigenvalueDecomposition = new EigenvalueDecomposition(new Matrix(Aa));
//        roots = ComplexArrays.zip(eigenvalueDecomposition.getRealEigenvalues(), eigenvalueDecomposition.getImagEigenvalues());
//        coefficients = new Polynomial(roots).getCoefficients();
//        System.out.println(Arrays.toString(coefficients));
//
//        System.out.println(Arrays.toString(new Polynomial(new Complex[]{Complex.fromReal(-1)}).getCoefficients()));
//
//        TransferFunction tf = new TransferFunction(new double[]{1}, new double[]{1});
//        System.out.println(tf);
//
//        StateSpace ss = tf.toStateSpace();
//        System.out.println(ss);

        Test2();
    }

    public static void Test2() {
        double value = -5.35;
        double value1 = 0.0;
        double value2 = -0.0;
        double value3 = Double.NaN;
        double value4 = Double.NEGATIVE_INFINITY;
        double value5 = Double.POSITIVE_INFINITY;

        MathETK.FRexpResult frexp = frexp(value);
        System.out.println("Mantissa: " + frexp.mantissa);
        System.out.println("Exponent: " + frexp.exponent);
        System.out.println("Original value was: " + value);
        System.out.println(frexp.mantissa + " * 2^" + frexp.exponent + " = ");
        System.out.println(frexp.mantissa * (1 << frexp.exponent));
        System.out.println();

        frexp = frexp(value1);
        System.out.println("Mantissa: " + frexp.mantissa);
        System.out.println("Exponent: " + frexp.exponent);
        System.out.println("Original value was: " + value1);
        System.out.println(frexp.mantissa + " * 2^" + frexp.exponent + " = ");
        System.out.println(frexp.mantissa * (1 << frexp.exponent));
        System.out.println();

        frexp = frexp(value2);
        System.out.println("Mantissa: " + frexp.mantissa);
        System.out.println("Exponent: " + frexp.exponent);
        System.out.println("Original value was: " + value2);
        System.out.println(frexp.mantissa + " * 2^" + frexp.exponent + " = ");
        System.out.println(frexp.mantissa * (1 << frexp.exponent));
        System.out.println();


        frexp = frexp(value3);
        System.out.println("Mantissa: " + frexp.mantissa);
        System.out.println("Exponent: " + frexp.exponent);
        System.out.println("Original value was: " + value3);
        System.out.println(frexp.mantissa + " * 2^" + frexp.exponent + " = ");
        System.out.println(frexp.mantissa * (1 << frexp.exponent));
        System.out.println();

        frexp = frexp(value4);
        System.out.println("Mantissa: " + frexp.mantissa);
        System.out.println("Exponent: " + frexp.exponent);
        System.out.println("Original value was: " + value4);
        System.out.println(frexp.mantissa + " * 2^" + frexp.exponent + " = ");
        System.out.println(frexp.mantissa * (1 << frexp.exponent));
        System.out.println();

        frexp = frexp(value5);
        System.out.println("Mantissa: " + frexp.mantissa);
        System.out.println("Exponent: " + frexp.exponent);
        System.out.println("Original value was: " + value5);
        System.out.println(frexp.mantissa + " * 2^" + frexp.exponent + " = ");
        System.out.println(frexp.mantissa * (1 << frexp.exponent));
    }

    /**
     * Characteristic polynomial of matrix.
     *
     * @return The Characteristic polynomial of the Matrix.
     */
    public double[] poly() {
        EigenvalueDecomposition eig = this.eig();
        Complex[] roots = ComplexArrays.zip(eig.getRealEigenvalues(), eig.getImagEigenvalues());
        return new Polynomial(roots).getCoefficients();
    }

}
