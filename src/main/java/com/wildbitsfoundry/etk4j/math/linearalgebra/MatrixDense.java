package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.Arrays;
import java.util.Random;

import static com.wildbitsfoundry.etk4j.math.MathETK.frExp;

public class MatrixDense extends Matrix {
    private double[] data;

    /**
     * Constructs a {@code Matrix}.
     *
     * @param rows The number of rows of the {@code Matrix}..
     * @param cols The number of columns of the {@code Matrix}..
     */
    public MatrixDense(int rows, int cols) {
        super(rows, cols);
        this.data = new double[rows * cols];
    }

    /***
     * Construct {@code Matrix}  from column packed data. Column packed data can be represented as:
     * <pre>
     *     For a given Matrix:
     *     A = 1 2 3
     *         4 5 6
     *         7 8 9
     *     The column packed data will be [1, 4, 7, 2, 5, 8, 3, 6, 9]
     * </pre>
     * No checks are done to ensure that the data is not null so proceed with caution.
     * @param data The column packed data.
     * @param rows The number of rows of the {@code Matrix}.
     */
    public MatrixDense(double[] data, int rows) {
        this.rows = rows;
        cols = (this.rows != 0 ? data.length / this.rows : 0);
        if (this.rows * cols != data.length) {
            throw new IllegalArgumentException("Array length must be a multiple of rows.");
        }

        int dim = this.rows * cols;
        this.data = new double[dim];
        for (int i = 0; i < this.rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                this.data[i * cols + j] = data[i + j * rows];
            }
        }
    }

    /**
     * Construct matrix from a non-jagged 2D array of values.
     *
     * @param data The 2D array of values used to populate the {@code Matrix} internal storage.
     */
    public MatrixDense(double[][] data) {
        if(data == null) {
            rows = 0;
            cols = 0;
            this.data = null;
        } else if(data.length == 0) {
            rows = 0;
            cols = 0;
            this.data = new double[0];
        } else {
            rows = data.length;
            cols = data[0].length;
            this.data = DoubleArrays.flatten(data);
        }
    }

    /***
     * Construct {@code Matrix}  from row packed data. Row packed data can be represented as:
     * <pre>
     *     For a given Matrix:
     *     A = 1 2 3
     *         4 5 6
     *         7 8 9
     *     The row packed data will be [1, 2, 3, 4, 5, 6, 7, 8, 9]
     * </pre>
     * No checks are done to ensure that the data is not null so proceed with caution.
     * @param data The row packed data.
     * @param rows The number of rows of the {@code Matrix}.
     * @param cols The number of columns of the {@code Matrix}.
     */
    public MatrixDense(double[] data, int rows, int cols) {
        this.rows = rows;
        this.cols = cols;
        this.data = data;
    }

    /**
     * Copy constructor.
     *
     * @param matrix The argument {@code Matrix} to be copied.
     */
    public MatrixDense(MatrixDense matrix) {
        rows = matrix.rows;
        cols = matrix.cols;
        data = new double[rows * cols];
        System.arraycopy(matrix.data, 0, this.data, 0, this.rows * this.cols);
    }

    /**
     * Constructs a {@code Matrix} filled with a predefined value.
     *
     * @param rows The number of rows of the {@code Matrix}.
     * @param cols The number of columns of the {@code Matrix}.
     * @param val  The argument used to fill up the {@code Matrix}.
     */
    public MatrixDense(int rows, int cols, double val) {
        this.rows = rows;
        this.cols = cols;
        data = new double[this.rows * this.cols];
        Arrays.fill(data, val);
    }

    /***
     * Deep copy.
     * @return A deep copy of the {@code Matrix}.
     */
    public MatrixDense copy() {
        double[] data = Arrays.copyOf(this.data, this.data.length);
        return new MatrixDense(data, this.rows, this.cols);
    }

    /***
     * Creates a new {@link MatrixDense} from a 2D array of values
     * @param A the input array
     * @return a new {@link MatrixDense} created from the input array
     */
    public static MatrixDense from2DArray(double[][] A) {
        return new MatrixDense(A);
    }

    /***
     * Get sub-matrix.
     *
     * @param row0 The initial row index.
     * @param row1 The final row index.
     * @param col0 The initial column index.
     * @param col1 The final column index.
     * @return {@code A(row0 : row1, col0 : col1)}.
     */
    public MatrixDense subMatrix(int row0, int row1, int col0, int col1) {
        if (row0 < 0 || row1 < 0) {
            throw new IllegalArgumentException("The row indexes row0 and row1 must be greater than zero.");
        }
        if (col0 < 0 || col1 < 0) {
            throw new IllegalArgumentException("The column indexes col0 and col1 must be greater than zero.");
        }
        if (row1 >= rows) {
            throw new ArrayIndexOutOfBoundsException("The final row index cannot be greater than the number of rows in the Matrix.");
        }
        if (col1 >= cols) {
            throw new ArrayIndexOutOfBoundsException("The final column index cannot be greater than the number of columns in the Matrix.");
        }
        int rowDim = row1 - row0 + 1;
        int colDim = col1 - col0 + 1;
        double[] data = new double[rowDim * colDim];
        for (int i = row0; i <= row1; ++i) {
            for (int j = col0; j <= col1; ++j) {
                data[(i - row0) * colDim + (j - col0)] = this.data[i * cols + j];
            }
        }
        return new MatrixDense(data, rowDim, colDim);
    }

    /***
     * Get sub-matrix
     *
     * @param rows The array of row indices.
     * @param col0 The initial column index.
     * @param col1 The final column index.
     * @return {@code A(rows(:), col0 : col1)}.
     */
    public MatrixDense subMatrix(int[] rows, int col0, int col1) {
        if (col0 < 0 || col1 < 0) {
            throw new IllegalArgumentException("The column indexes col0 and col1 must be greater than zero.");
        }
        if (Arrays.stream(rows).anyMatch(i -> i >= this.rows || i < 0)) {
            throw new ArrayIndexOutOfBoundsException("The row indexes cannot be greater than the number of rows in the Matrix. and must be greater than zero.");
        }
        if (col1 >= cols) {
            throw new ArrayIndexOutOfBoundsException("The final column index cannot be greater than the number of columns in the Matrix.");
        }
        int rowDim = rows.length;
        int colDim = col1 - col0 + 1;
        double[] data = new double[rowDim * colDim];
        for (int i = 0; i < rowDim; ++i) {
            for (int j = 0; j < colDim; ++j) {
                data[i * colDim + j] = this.data[rows[i] * cols + (j + col0)];
            }
        }
        return new MatrixDense(data, rowDim, colDim);
    }

    /***
     * Get sub-matrix.
     *
     * @param row0 The initial row index.
     * @param row1 The final row index.
     * @param cols The array of column indices.
     * @return {@code A(row0 : row1, cols (:))}.
     */
    public MatrixDense subMatrix(int row0, int row1, int[] cols) {
        if (row0 < 0 || row1 < 0) {
            throw new ArrayIndexOutOfBoundsException("The row indexes row0 and row1 must be greater than zero.");
        }
        if (Arrays.stream(cols).anyMatch(i -> i >= this.cols || i < 0)) {
            throw new ArrayIndexOutOfBoundsException("The column indexes cannot be greater than the number of columns in the Matrix and must be greater than zero.");
        }
        if (row1 >= rows) {
            throw new ArrayIndexOutOfBoundsException("The final ros index cannot be greater than the number of rows in the Matrix.");
        }
        int rowDim = row1 - row0 + 1;
        int colDim = cols.length;
        double[] data = new double[rowDim * colDim];
        for (int i = 0; i < rowDim; ++i) {
            for (int j = 0; j < colDim; ++j) {
                data[i * colDim + j] = this.data[(i + row0) * this.cols + cols[j]];
            }
        }
        return new MatrixDense(data, rowDim, colDim);
    }

    /***
     * Get sub-matrix.
     *
     * @param rows The array or row indices.
     * @param cols The array of column indices
     * @return {@code A(rows(:), cols(:))}.
     */
    public MatrixDense subMatrix(int[] rows, int[] cols) {
        if (Arrays.stream(rows).anyMatch(i -> i >= this.rows || i < 0)) {
            throw new ArrayIndexOutOfBoundsException("The row indexes cannot be greater than the number of rows in the Matrix. and must be greater than zero.");
        }
        if (Arrays.stream(cols).anyMatch(i -> i >= this.cols || i < 0)) {
            throw new ArrayIndexOutOfBoundsException("The column indexes cannot be greater than the number of columns in the Matrix and must be greater than zero.");
        }
        int rowDim = rows.length;
        int colDim = cols.length;
        double[] data = new double[rowDim * colDim];
        for (int i = 0; i < rowDim; ++i) {
            for (int j = 0; j < colDim; ++j) {
                data[i * colDim + j] = this.data[rows[i] * this.cols + cols[j]];
            }
        }
        return new MatrixDense(data, rowDim, colDim);
    }

    /**
     * Retrieve value from {@code Matrix} at a given position.
     *
     * @param row The row index.
     * @param col The column index.
     * @return The value at {@code A(i, j)}.
     */
    public double get(int row, int col) {
        if (row < 0 || row >= rows || col < 0 || col >= cols)
            throw new ArrayIndexOutOfBoundsException("Outside of matrix bounds");
        return data[row * cols + col];
    }

    @Override
    public double unsafeGet(int row, int col) {
        return data[row * cols + col];
    }

    /**
     * Set the value of the {@code Matrix} at a given position.
     *
     * @param row   The row index.
     * @param col   The column index.
     * @param val The value used to set {@code A(i, j) = val}.
     */
    public void set(int row, int col, double val) {
        if (row < 0 || row >= rows || col < 0 || col >= cols)
            throw new ArrayIndexOutOfBoundsException("Outside of matrix bounds");
        data[row * cols + col] = val;
    }

    @Override
    public void unsafeSet(int row, int col, double val) {
        data[row * cols + col] = val;
    }

    /**
     * {@code Matrix} determinant.
     *
     * @return The determinant of the squared {@code Matrix}.
     */
    public double det() {
        return new LUDecompositionDense(this).det();
    }

    /**
     * Set a sub-matrix.
     *
     * @param i0 Initial row index
     * @param i1 Final row index
     * @param j0 Initial column index
     * @param j1 Final column index
     * @param X  A(i0:i1,j0:j1)
     */
    public void setMatrix(int i0, int i1, int j0, int j1, MatrixDense X) {
        if(i0 < 0 || i0 >= rows) {
            throw new ArrayIndexOutOfBoundsException("The initial row index i0 must be greater than zero and less than the number of rows.");
        }
        if(i1 < 0 || i1 >= rows) {
            throw new ArrayIndexOutOfBoundsException("The final row index i1 must be greater than zero and less than the number of rows.");
        }
        if(j0 < 0 || j0 >= cols) {
            throw new ArrayIndexOutOfBoundsException("The initial column index j0 must be greater than zero and less than the number of columns.");
        }
        if(j1 < 0 || j1 >= cols) {
            throw new ArrayIndexOutOfBoundsException("The final column index j1 must be greater than zero and less than the number of columns.");
        }
        for (int i = i0; i <= i1; i++) {
            for (int j = j0; j <= j1; j++) {
                data[i * cols + j] = X.get(i - i0, j - j0);
            }
        }
    }

    /**
     * Set a sub-matrix.
     *
     * @param r Array of row indices.
     * @param c Array of column indices.
     * @param X {@code A(r(:),c(:))}
     */
    public void setMatrix(int[] r, int[] c, MatrixDense X) {
        if(Arrays.stream(r).anyMatch(i -> i < 0 || i >= rows)) {
            throw new ArrayIndexOutOfBoundsException("The row indexes must be greater than zero and less than the number of rows.");
        }
        if(Arrays.stream(c).anyMatch(j -> j < 0 || j >= cols)) {
            throw new ArrayIndexOutOfBoundsException("The column indexes must be greater than zero and less than the number of columns.");
        }
        for (int i = 0; i < r.length; i++) {
            for (int j = 0; j < c.length; j++) {
                data[r[i] * cols + c[j]] = X.get(i, j);
            }
        }
    }

    /**
     * Set a sub-matrix.
     *
     * @param r  Array of row indices.
     * @param j0 Initial column index
     * @param j1 Final column index
     * @param X  A(r(:),j0:j1)
     */
    public void setMatrix(int[] r, int j0, int j1, MatrixDense X) {
        if(Arrays.stream(r).anyMatch(i -> i < 0 || i >= rows)) {
            throw new ArrayIndexOutOfBoundsException("The row indexes must be greater than zero and less than the number of rows.");
        }
        if(j0 < 0 || j0 >= cols) {
            throw new ArrayIndexOutOfBoundsException("The initial column index j0 must be greater than zero and less than the number of columns.");
        }
        if(j1 < 0 || j1 >= cols) {
            throw new ArrayIndexOutOfBoundsException("The final column index j1 must be greater than zero and less than the number of columns.");
        }
        for (int i = 0; i < r.length; i++) {
            for (int j = j0; j <= j1; j++) {
                data[r[i] * cols + j] = X.get(i, j - j0);
            }
        }
    }

    /**
     * Set a sub-matrix.
     *
     * @param i0 Initial row index
     * @param i1 Final row index
     * @param c  Array of column indices.
     * @param X  A(i0:i1,c(:))
     */
    public void setMatrix(int i0, int i1, int[] c, MatrixDense X) {
        if(i0 < 0 || i0 >= rows) {
            throw new ArrayIndexOutOfBoundsException("The initial row index i0 must be greater than zero and less than the number of rows.");
        }
        if(i1 < 0 || i1 >= rows) {
            throw new ArrayIndexOutOfBoundsException("The final row index i1 must be greater than zero and less than the number of rows.");
        }
        if(Arrays.stream(c).anyMatch(j -> j < 0 || j >= cols)) {
            throw new ArrayIndexOutOfBoundsException("The column indexes must be greater than zero and less than the number of columns.");
        }
        for (int i = i0; i <= i1; i++) {
            for (int j = 0; j < c.length; j++) {
                data[i * cols + c[j]] = X.get(i - i0, j);
            }
        }
    }

    /**
     * Matrix rank.
     *
     * @return The effective numerical rank, obtained from SVD.
     */
    public int rank() {
        return new SingularValueDecompositionDense(this).rank();
    }

    /**
     * Matrix condition (2 norm)
     *
     * @return ratio of largest to the smallest singular value.
     */

    public double cond() {
        return new SingularValueDecompositionDense(this).cond();
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

    /**
     * Cofactor {@code Matrix}
     *
     * @return The cofactor {@code Matrix}.
     * @see <a href="https://mathworld.wolfram.com/Cofactor.html">Cofactor Matrix.</a>
     */
    public MatrixDense cofactor() {
        int dim = this.rows;
        if (!this.isSquare()) {
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
        return new MatrixDense(cofactor);
    }

    /**
     * {@code Matrix} minor.
     *
     * @param row The row index.
     * @param col The column index.
     * @return The minor of the {@code Matrix} by omitting the {@code row} and {@code column} indexes.
     * @see <a href="https://mathworld.wolfram.com/Minor.html">Matrix minor</a>
     */
    public MatrixDense minor(int row, int col) {
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
        return new MatrixDense(minor);
    }

    /**
     * Adjoint (adjugate) {@code Matrix}.
     *
     * @return The {@link MatrixDense#transpose()} of the {@link MatrixDense#adjoint()} {@code Matrix}.
     * @see <a href="https://en.wikipedia.org/wiki/Adjugate_matrix">Adjoint matrix</a>
     */
    public MatrixDense adjoint() {
        return this.cofactor().transpose();
    }

    // region decompositions

    /**
     * LU decomposition of the {@code Matrix}.
     *
     * @return The {@link LUDecompositionDense} of the {@code Matrix}.
     * @see <a href="https://en.wikipedia.org/wiki/LU_decomposition">LU Decomposition</a>
     */
    public LUDecompositionDense LU() {
        return new LUDecompositionDense(this);
    }

    /**
     * QR decomposition of the {@code Matrix}.
     *
     * @return The {@link QRDecompositionDense} of the {@code Matrix}.
     * @see <a href="https://en.wikipedia.org/wiki/QR_decomposition">QR Decomposition</a>
     */
    public QRDecompositionDense QR() {
        return new QRDecompositionDense(this);
    }

    /**
     * Cholesky decomposition of the {@code Matrix}.
     *
     * @return The {@link CholeskyDecompositionDense} of the {@code Matrix}.
     * @see <a href="https://en.wikipedia.org/wiki/Cholesky_decomposition">Cholesky Decomposition</a>
     */
    public CholeskyDecompositionDense Chol() {
        return new CholeskyDecompositionDense(this);
    }

    /**
     * Singular Value Decomposition of the {@code Matrix}.
     *
     * @return The {@link SingularValueDecompositionDense} of the {@code Matrix}.
     * @see <a href="https://en.wikipedia.org/wiki/Singular_value_decomposition">Singular Value Decomposition</a>
     */
    public SingularValueDecompositionDense SVD() {
        return new SingularValueDecompositionDense(this);
    }
    // endregion

    /***
     * Get matrix diagonal
     *
     * @return Array containing the diagonal of the matrix
     */
    public double[] diag() {
        if (!this.isSquare()) {
            throw new RuntimeException("Matrix is not squared");
        }
        final int dim = rows;
        double[] diag = new double[dim];
        for (int i = 0; i < dim; ++i) {
            diag[i] = data[i * dim + i];
        }
        return diag;
    }

    /**
     * Inverse of the {@code Matrix}.
     *
     * @return {@code A<sup>-1</sup>}.
     */
    public MatrixDense inv() {
        return this.solve(MatrixDense.Factory.identity(rows));
    }

    /**
     * Is the {@code Matrix} empty.
     *
     * @return {@code true} if the number of columns and rows are equal to zero or the internal data is null or the
     * length of the internal data is equal to zero. Returns {@code false} otherwise.
     */
    public boolean isEmpty() {
        if ((rows == 0 && cols == 0) || data == null || data.length == 0) {
            return true;
        }
        return false;
    }

    /**
     * {@code Matrix} transpose.
     *
     * @return {@code A<sup>T</sup>}.
     * @see <a href="https://en.wikipedia.org/wiki/Transpose">Matrix transpose</a>
     */
    public MatrixDense transpose() {
        if (this.isEmpty()) {
            return MatrixDense.Factory.empty();
        }
        double[] result = new double[rows * cols];
        final int trows = cols;
        final int tcols = rows;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result[j * tcols + i] = data[i * cols + j];
            }
        }
        return new MatrixDense(result, trows, tcols);
    }

    // region norms

    /***
     * One norm
     *
     * @return maximum column sum.
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
        return new SingularValueDecompositionDense(this).norm2();
    }

    /***
     * Infinity norm
     *
     * @return maximum row sum.
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
     * @return square root of the sum of squares of all elements.
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
    public MatrixDense uminus() {
        final int m = rows;
        final int n = cols;
        final int length = m * n;
        double[] x = new double[length];
        for (int i = 0; i < length; ++i) {
            x[i] = -data[i];
        }
        return new MatrixDense(x, m, n);
    }

    // region arithmetic operations

    /**
     * Matrix addition.
     * @param M The {@code Matrix} to add.
     * @return {@code A + M}.
     */
    public MatrixDense add(MatrixDense M) {
        checkMatrixDimensions(M);
        double[] result = new double[this.rows * this.cols];
        for (int i = 0; i < this.rows * this.cols; ++i) {
            result[i] = this.data[i] + M.data[i];
        }
        return new MatrixDense(result, rows, cols);
    }

    /**
     * Matrix addition in place. This is equivalent to {@code A += B}.
     * @param B The {@code Matrix} to add.
     */
    public void addEquals(MatrixDense B) {
        checkMatrixDimensions(B);
        final int length = rows * cols;
        for (int i = 0; i < length; ++i) {
            data[i] += B.data[i];
        }
    }

    /**
     * Matrix subtraction.
     * @param M The {@code Matrix} to subtract.
     * @return {@code A - M}.
     */
    public MatrixDense subtract(MatrixDense M) {
        checkMatrixDimensions(M);
        double[] result = new double[this.rows * this.cols];
        for (int i = 0; i < this.rows * this.cols; ++i) {
            result[i] = this.data[i] - M.data[i];
        }
        return new MatrixDense(result, rows, cols);
    }

    /**
     * Matrix subtraction in place. This is equivalent to {@code A -= B}.
     * @param B The {@code Matrix} to subtract.
     */
    public void subtractEquals(MatrixDense B) {
        checkMatrixDimensions(B);
        final int length = rows * cols;
        for (int i = 0; i < length; ++i) {
            data[i] -= B.data[i];
        }
    }

    /**
     * Element-by-element multiplication, C = A.*B
     *
     * @param B another matrix
     * @return A.*B
     */
    public MatrixDense arrayMultiply(MatrixDense B) {
        checkMatrixDimensions(B);
        return new MatrixDense(DoubleArrays.multiplyElementWise(data, B.data), rows, cols);
    }

    /**
     * Element-by-element multiplication in place, A = A.*B
     *
     * @param B another matrix
     */
    public void arrayMultiplyEquals(MatrixDense B) {
        checkMatrixDimensions(B);
        DoubleArrays.multiplyElementWiseInPlace(data, B.data);
    }

    /**
     * Element-by-element right division, C = A./B
     *
     * @param B another matrix
     * @return A./B
     */
    public MatrixDense arrayRightDivide(MatrixDense B) {
        checkMatrixDimensions(B);
        return new MatrixDense(DoubleArrays.divideElementWise(data, B.data), rows, cols);
    }

    /**
     * Element-by-element right division in place, A = A./B
     *
     * @param B another matrix
     */
    public void arrayRightDivideEquals(MatrixDense B) {
        checkMatrixDimensions(B);
        DoubleArrays.divideElementWiseInPlace(data, B.data);
    }

    /**
     * Element-by-element left division, C = A.\B
     *
     * @param B another matrix
     * @return A.\B
     */
    public MatrixDense arrayLeftDivide(MatrixDense B) {
        checkMatrixDimensions(B);
        if (data.length != B.data.length) {
            throw new IllegalArgumentException("a and b dimensions must match");
        }

        final int length = data.length;
        double[] data = Arrays.copyOf(this.data, length);
        for (int i = 0; i < length; ++i) {
            data[i] = B.data[i] / data[i];
        }
        return new MatrixDense(data, rows, cols);
    }

    /**
     * Element-by-element left division in place, A = A.\B
     *
     * @param B another matrix
     */
    public void arrayLeftDivideEquals(MatrixDense B) {
        checkMatrixDimensions(B);
        if (data.length != B.data.length) {
            throw new IllegalArgumentException("a and b dimensions must match");
        }

        final int length = data.length;
        for (int i = 0; i < length; ++i) {
            data[i] = B.data[i] / data[i];
        }
    }

    /**
     * Multiply a matrix by a scalar, C = s*A
     *
     * @param s scalar
     * @return s*A
     */
    public MatrixDense multiply(double s) {
        return new MatrixDense(DoubleArrays.multiplyElementWise(data, s), rows, cols);
    }

    /**
     * Multiply a matrix by a complex scalar, C = s*A
     *
     * @param s scalar
     * @return s*A
     */
    public ComplexMatrixDense multiply(Complex s) {
        return new ComplexMatrixDense(ComplexArrays.multiplyElementWise(data, s), rows, cols);
    }

    /**
     * Multiply a matrix by a scalar in place, A = s*A
     *
     * @param s scalar
     */
    public void multiplyEquals(double s) {
        DoubleArrays.multiplyElementWiseInPlace(data, s);
    }

    /**
     * Matrix multiplication.
     * @param B The matrix to multiply.
     * @return {@code A * B}
     */
    public MatrixDense multiply(MatrixDense B) {
        MatrixDense c = new MatrixDense(0, 0);
        multiplyOp(this, B, c);
        return c;
    }

    /**
     * Matrix multiplication in place. This is equivalent to {@code A *= B}
     * @param B The {@code Matrix} to multiply.
     */
    public void multiplyEquals(MatrixDense B) {
        multiplyOp(this, B, this);
    }

    private static void multiplyOp(MatrixDense a, MatrixDense b, MatrixDense c) {
        if (b.rows != a.cols) {
            throw new IllegalArgumentException("Matrix inner dimensions must agree. Check that the number of " +
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

    /**
     * Complex {@code Matrix} multiplication.
     * @param matrix The complex {@code Matrix} to multiply;
     * @return {@code A * B}.
     */
    public ComplexMatrixDense multiply(ComplexMatrixDense matrix) {
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
        return new ComplexMatrixDense(result, rows, bCols);
    }

    // endregion

    /**
     * Solve system of linear equations. Three different algorithms are used depending on the shape of the matrix:
     * <pre>
     *     LU Decomposition if the matrix is squared.
     *     QR if the matrix is thin in other words it has more rows than columns. (Overdetermined system)
     *     Transpose QR if the matrix is short and wide in other words it has more columns than rows. (Under-determined system)
     * </pre>
     *
     * @param b The solution {@link Matrix}.
     * @return The solution to {@code Ax = b}
     */
    public MatrixDense solve(MatrixDense b) {
        if (rows == cols) { // Matrix is Squared
            return new LUDecompositionDense(this).solve(b);
        } else if (rows > cols) { // Matrix is tall and narrow (Overdetermined system)
            return new QRDecompositionDense(this).solve(b);
        } else { // Matrix is short and wide (Under-determined system)
            // Could use QR for matrices that are not rank deficient. Let's go
            // with pinv since we don't know what the input matrix looks like
            return this.pinv().multiply(b);
        }
    }

    // TODO document and write test
    public MatrixDense solve(double[] b) {
        return solve(new MatrixDense(b, b.length));
    }

    /**
     * Solve X*A = B, which is also A'*X' = B'
     *
     * @param B right hand side
     * @return solution if A is square, least squares solution otherwise.
     */
    public MatrixDense transposeSolve(MatrixDense B) {
        return transpose().solve(B.transpose());
    }

    public void appendRows(int count) {
        rows += count;
        final int newSize = rows * cols;
        data = Arrays.copyOf(data, newSize);
    }

    /**
     * Pseudo inverse of the {@code Matrix}.
     *
     * @return {@code A<sup>+</sup>}.
     */
    public MatrixDense pinv() {
        int rows = this.rows;
        int cols = this.cols;

        if (rows < cols) {
            MatrixDense result = this.transpose().pinv();
            if (result != null) {
                result = result.transpose();
            }
            return result;
        }

        SingularValueDecompositionDense svdX = this.SVD();
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
        MatrixDense U = svdX.getU();
        MatrixDense V = svdX.getV();
        int min = Math.min(cols, U.cols);
        double[][] inverse = new double[cols][rows];
        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < U.rows; j++) {
                for (int k = 0; k < min; k++) {
                    inverse[i][j] += V.get(i, k) * singularValueReciprocals[k] * U.get(j, k);
                }
            }
        }
        return new MatrixDense(inverse);
    }

    @Override
    public String toString() {
        if(isEmpty()) {
            return "[]";
        }
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

    /**
     * Eigenvalue decomposition. The {@link Matrix} is balanced ({@link MatrixDense#balance()}) prior to the decomposition.
     *
     * @return The {@link EigenvalueDecompositionDense} of the {@code Matrix}.
     */
    public EigenvalueDecompositionDense eig() {
        return new EigenvalueDecompositionDense(this);
    }

    /**
     * Eigenvalue decomposition with optional pre balance.
     *
     * @param balance If this flag is set to {@code true}, the {@code Matrix} is balanced prior to the decomposition.
     * @return The {@link EigenvalueDecompositionDense} of the {@code Matrix}.
     */
    public EigenvalueDecompositionDense eig(boolean balance) {
        return new EigenvalueDecompositionDense(this, balance);
    }

    /**
     * Balances the matrix using the algorithm by Parlett and Reinsch with norm -1.
     * References:
     * <pre>
     * <a href="http://www.netlib.org/eispack/balanc.f">http://www.netlib.org/eispack/balanc.f</a>
     * <a href="https://arxiv.org/pdf/1401.5766.pdf">https://arxiv.org/pdf/1401.5766.pdf - Algorithm #2</a>
     * </pre>
     *
     * @return A balanced copy of the {@code Matrix}.
     * @see <a href="https://stackoverflow.com/a/43169781/6383857">Balancing a matrix</a>
     * @see <a href="https://www.mathworks.com/help/matlab/ref/balance.html">Matrix balance</a>
     */
    public MatrixDense balance() {
        if (!this.isSquare()) {
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
        return new MatrixDense(data, rows, rows);
    }

    /**
     * Number of rows.
     *
     * @return The number of rows of the {@code Matrix}.
     */
    public int getRowCount() {
        return rows;
    }

    /**
     * Number of columns.
     *
     * @return The number of columns of the {@code Matrix}.
     */
    public int getColumnCount() {
        return cols;
    }

    /**
     * Retrieve a row.
     *
     * @param row The index of the column to retrieve.
     * @return The values of the row at the specified index.
     */
    public double[] getRow(int row) {
        double[] result = new double[cols];
        int rowIndex = row * cols;
        System.arraycopy(data, rowIndex + 0, result, 0, cols);
        return result;
    }

    /**
     * Retrieve a column.
     *
     * @param col The index of the column to retrieve.
     * @return The values of the column at the specified index.
     */
    public double[] getCol(int col) {
        if (col < 0 || col >= cols) {
            throw new IllegalArgumentException("The column index col must be greater than zero and less than the number of columns.");
        }
        double[] result = new double[rows];
        for (int i = 0; i < rows; ++i) {
            result[i] = data[i * rows + col];
        }
        return result;
    }

    /**
     * Retrieves the copy of the row packed data of the {@code Matrix}.
     * <pre>
     *     For a given Matrix:
     *     A = 1 2 3
     *         4 5 6
     *         7 8 9
     *     The row packed data will be [1, 2, 3, 4, 5, 6, 7, 8, 9]
     * </pre>
     *
     * @return A copy of the internal data of the {@code Matrix}.
     */
    public double[] getArrayCopy() {
        return Arrays.copyOf(data, data.length);
    }

    /**
     * Retrieves the row packed data of the {@code Matrix}. Changes to this array will change the contents of the {@code Matrix}.
     * <pre>
     *     For a given Matrix:
     *     A = 1 2 3
     *         4 5 6
     *         7 8 9
     *     The row packed data will be [1, 2, 3, 4, 5, 6, 7, 8, 9]
     * </pre>
     *
     * @return The internal data of the {@code Matrix}.
     */
    public double[] getArray() {
        return data;
    }

    /**
     * Retrieve {@link Matrix} data as 2D array.
     *
     * @return A 2D array copy of the internal data of the Matrix.
     */
    public double[][] getAs2DArray() {
        double[][] data = new double[rows][cols];
        for (int i = 0; i < rows; ++i) {
            data[i] = Arrays.copyOfRange(this.data, i * cols, i * cols + cols);
        }
        return data;
    }

    /**
     * Make a one-dimensional column packed copy of the internal array.
     * <pre>
     *     For a given Matrix:
     *     A = 1 2 3
     *         4 5 6
     *         7 8 9
     *     The column packed data will be [1, 4, 7, 2, 5, 8, 3, 6, 9]
     * </pre>
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
        if (!(obj instanceof MatrixDense)) {
            return false;
        }
        MatrixDense other = (MatrixDense) obj;
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
     * <pre>
     *     For a given Matrix:
     *     A = 1 2 3
     *         4 5 6
     *         7 8 9
     *     The row packed data will be [1, 2, 3, 4, 5, 6, 7, 8, 9]
     * </pre>
     * @return Matrix elements packed in a one-dimensional array by rows.
     */
    public double[] getRowPackedCopy() {
        return this.getArrayCopy();
    }

    /**
     * Check if size(A) == size(B)
     **/
    void checkMatrixDimensions(MatrixDense B) {
        if (B.rows != rows || B.cols != cols) {
            throw new IllegalArgumentException("Matrix dimensions must agree.");
        }
    }

    /**
     * Set a row to a predefined set of values.
     *
     * @param i   The index of the row to set.
     * @param row The values to set.
     */
    public void setRow(int i, double[] row) {
        if (i < 0 || i >= rows) {
            throw new IndexOutOfBoundsException("The row index i is out of bounds, it must be greater than zero and less than the number of rows.");
        }
        if (row.length != cols) {
            throw new IllegalArgumentException("The length of the row cannot be greater than the number of columns.");
        }
        System.arraycopy(row, 0, data, i * cols, cols);
    }

    /**
     * Power of Matrix.
     *
     * @param n The power to raise the matrix to.
     * @return {@code A<sup>n</n>}.
     */
    public MatrixDense pow(int n) {
        if (!this.isSquare()) {
            throw new IllegalArgumentException("Matrix must be a square matrix.");
        }
        n = Math.abs(n);
        if (n == 0) {
            return MatrixDense.Factory.identity(rows);
        }
        MatrixDense a = n < 0 ? inv() : new MatrixDense(this);

        if (n == 1) {
            return a;
        } else if (n == 2) {
            a.multiplyEquals(a);
            return a;
        } else if (n == 3) {
            a.multiplyEquals(a.multiply(a));
            return a;
        }
        MatrixDense z = null, result = null;
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
     * Power of Matrix using the Eigenvalues method.
     *
     * @param n The power to raise the matrix to.
     * @return {@code A<sup>n</n>}.
     */
    public MatrixDense pow(double n) {
        if (!this.isSquare()) {
            throw new IllegalArgumentException("Matrix must be squared");
        }
        n = Math.abs(n);
        if (n == 0.0 || n == 1.0 || n == 2.0 || n == 3.0) {
            return pow((int) n);
        }
        MatrixDense a = n < 0 ? inv() : new MatrixDense(this);
        EigenvalueDecompositionDense eig = a.eig();
        MatrixDense D = eig.getD();
        for (int i = 0, j = 0; i < D.rows; ++i, ++j) {
            D.set(i, j, Math.pow(D.get(i, j), n));
        }
        MatrixDense result = eig.getV().multiply(D);
        result.multiplyEquals(eig.getV().inv());
        return result;
    }

    /**
     * {@code Matrix} exponential.
     *
     * @return {@code e<sup>A</sup>}.
     */
    public MatrixDense expm() {
        MathETK.FRexpResult result = frExp(this.normInf());
        double s = Math.max(0.0, result.exponent + 1);
        MatrixDense A = this.multiply(1 / Math.pow(2.0, s));

        MatrixDense X = new MatrixDense(A);
        double c = 0.5;
        MatrixDense E = MatrixDense.Factory.identity(A.rows, A.cols).add(A.multiply(c));
        MatrixDense D = MatrixDense.Factory.identity(A.rows, A.cols).subtract(A.multiply(c));
        double q = 6.0;
        boolean p = true;
        for (int k = 2; k <= q; ++k) {
            c = c * (q - k + 1) / (k * (2 * q - k + 1));
            X.multiplyEquals(A);
            MatrixDense cX = X.multiply(c);
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

    /**
     * Characteristic polynomial of matrix.
     *
     * @return The Characteristic polynomial of the Matrix.
     */
    public double[] poly() {
        EigenvalueDecompositionDense eig = this.eig();
        Complex[] roots = ComplexArrays.zip(eig.getRealEigenvalues(), eig.getImagEigenvalues());
        return new Polynomial(roots).getCoefficients();
    }

    public static final class Factory {
        private Factory() {}

        /**
         * Vandermonde {@code Matrix}.
         * @param x The coefficient values.
         * @param rows The number of rows.
         * @param cols The number of columns.
         * @return {@code Vandermonde(rows, cols)}.
         * @see <a href="https://en.wikipedia.org/wiki/Vandermonde_matrix">Vandermonde matrix</a>
         */
        public static MatrixDense vandermonde(double[] x, int rows, int cols) {
            double[][] V = new double[rows][cols];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    V[i][j] = Math.pow(x[i], j);
                }
            }
            return new MatrixDense(V);
        }

        /**
         * Identity {@code Matrix}.
         * @param rows The number of rows.
         * @param cols The number of columns.
         * @return {@code identity(rows, cols)}.
         */
        public static MatrixDense identity(int rows, int cols) {
            double[] data = new double[rows * cols];
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    if (i == j) {
                        data[i * cols + j] = 1.0;
                    }
                }
            }
            return new MatrixDense(data, rows, cols);
        }

        /**
         * Identity {@code Matrix.}
         * @param n The number of rows and columns.
         * @return {@code identity(n, n)}.
         */
        public static MatrixDense identity(int n) {
            return MatrixDense.Factory.identity(n, n);
        }

        /**
         * Random {@code Matrix.}
         * @param n The number of rows and columns.
         * @return {@code random(n, n)}.
         */
        public static MatrixDense random(int n) {
            return random(n, n);
        }

        /**
         * Random {@code Matrix.}
         * @param rows The number of rows
         * @param cols The number of columns.
         * @return {@code random(rows, cols)}.
         */
        public static MatrixDense random(int rows, int cols) {
            Random rand = new Random();
            double[] data = new double[rows * cols];

            for (int i = 0; i < data.length; ++i) {
                data[i] = rand.nextDouble() * 100.0;
            }
            return new MatrixDense(data, rows, cols);
        }

        /**
         * Companion matrix.
         * @param coefficients The polynomial coefficients.
         * @param n The number of rows and columns.
         * @return The companion {@code Matrix} for the prescribed coefficients.
         * @see <a href="https://mathworld.wolfram.com/CompanionMatrix.html">Companion matrix</a>
         */
        public static MatrixDense companion(double[] coefficients, int n) {
            // Construct the companion matrix
            MatrixDense c = new MatrixDense(n, n);

            double a = 1.0 / coefficients[0];
            for (int i = 0; i < n; i++) {
                c.set(0, n - 1 - i, -coefficients[n - i] * a);
            }
            for (int i = 1; i < n; i++) {
                c.set(i, i - 1, 1);
            }
            return c;
        }

        /**
         * Create empty {@code Matrix}.
         * @return A {@code Matrix} with zero rows and zero columns.
         */
        public static MatrixDense empty() {
            return new MatrixDense(0, 0);
        }

        /**
         * Magic {@code Matrix}.
         * @param n The number of rows and columns.
         * @return A magic {@code Matrix} of dimensions {@code n}.
         * @see <a href="https://www.mathworks.com/help/matlab/ref/magic.html">Magic matrix</a>
         */
        public static MatrixDense magic(int n) {
            MatrixDense magicMatrix;
            if (n == 1) {
                magicMatrix = new MatrixDense(n, n);
                magicMatrix.set(0, 0, 1.0);
            } else if (n == 2) {
                return empty();
            }
            double[][] M = new double[n][n];

            // Odd order
            if ((n % 2) == 1) {
                int a = (n + 1) / 2;
                int b = (n + 1);
                for (int j = 0; j < n; j++) {
                    for (int i = 0; i < n; i++) {
                        M[i][j] = n * ((i + j + a) % n) + ((i + 2 * j + b) % n) + 1;
                    }
                }
                // Doubly Even Order
            } else if ((n % 4) == 0) {
                for (int j = 0; j < n; j++) {
                    for (int i = 0; i < n; i++) {
                        if (((i + 1) / 2) % 2 == ((j + 1) / 2) % 2) {
                            M[i][j] = n * n - n * i - j;
                        } else {
                            M[i][j] = n * i + j + 1;
                        }
                    }
                }
                // Singly Even Order
            } else {
                int p = n / 2;
                int k = (n - 2) / 4;
                MatrixDense A = magic(p);
                for (int j = 0; j < p; j++) {
                    for (int i = 0; i < p; i++) {
                        double aij = A.get(i, j);
                        M[i][j] = aij;
                        M[i][j + p] = aij + 2 * p * p;
                        M[i + p][j] = aij + 3 * p * p;
                        M[i + p][j + p] = aij + p * p;
                    }
                }
                for (int i = 0; i < p; i++) {
                    for (int j = 0; j < k; j++) {
                        double t = M[i][j];
                        M[i][j] = M[i + p][j];
                        M[i + p][j] = t;
                    }
                    for (int j = n - k + 1; j < n; j++) {
                        double t = M[i][j];
                        M[i][j] = M[i + p][j];
                        M[i + p][j] = t;
                    }
                }
                double t = M[k][0];
                M[k][0] = M[k + p][0];
                M[k + p][0] = t;

                t = M[k][k];
                M[k][k] = M[k + p][k];
                M[k + p][k] = t;
            }
            return new MatrixDense(M);
        }
    }
}