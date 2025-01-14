package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.Arrays;
import java.util.Random;
// TODO add balance method
public class ComplexMatrixDense extends ComplexMatrix {
    private Complex[] data;

    public ComplexMatrixDense(int rows, int cols) {
        super(rows, cols);
        this.data = new Complex[rows * cols];
    }
    // TODO add balance algo
    // TODO javadoc all file. Clean up the check in get set
    // add unsafe set and get
    public ComplexMatrixDense(Complex[] data, int rows) {
        this.rows = rows;
        cols = (this.rows != 0 ? data.length / this.rows : 0);
        if (this.rows * cols != data.length) {
            throw new IllegalArgumentException("Array length must be a multiple of rows.");
        }

        int dim = this.rows * cols;
        this.data = new Complex[dim];
        for (int i = 0; i < this.rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                this.data[i * cols + j] = data[i + j * rows];
            }
        }
    }

    public ComplexMatrixDense(Complex[][] data) {
        rows = data.length;
        cols = data[0].length;
        this.data = ComplexArrays.flatten(data);
    }

    // row packed
    // _rows = rows
    // _cols = cols
    // _data = data;
    public ComplexMatrixDense(Complex[] data, int rows, int cols) {
        this.rows = rows;
        this.cols = cols;
        this.data = data;
    }

    public ComplexMatrixDense(ComplexMatrixDense matrix) {
        rows = matrix.rows;
        cols = matrix.cols;
        data = new Complex[rows * cols];
        System.arraycopy(matrix.data, 0, this.data, 0, this.rows * this.cols);
    }

    public ComplexMatrixDense(int rows, int cols, Complex val) {
        this.rows = rows;
        this.cols = cols;
        data = new Complex[this.rows * this.cols];
        Arrays.fill(data, val);
    }

    public static ComplexMatrixDense from2DArray(Complex[][] A) {
        return new ComplexMatrixDense(A);
    }

    /***
     * Deep copy
     * @return A newly created {@link ComplexMatrixDense} with the values of this current one.
     */
    public ComplexMatrixDense copy() {
        Complex[] data = ComplexArrays.deepCopy(this.data);
        return new ComplexMatrixDense(data, this.rows, this.cols);
    }

    public static ComplexMatrixDense fromRealMatrix(MatrixDense m) {
        return new ComplexMatrixDense(ComplexArrays.fromReal(m.getArray()), m.getRowCount(), m.getColumnCount());
    }
    // region SubMatrix

    // endregion

    public Complex get(int i, int j) {
        // TODO clean up all these checks
        if (i < 0) {
            throw new ArrayIndexOutOfBoundsException("Index i cannot be less thant zero.");
        }
        if (i >= rows) {
            throw new ArrayIndexOutOfBoundsException(String.format("Index i: %d >= than number of rows: %d.", i, rows));
        }
        if (j < 0) {
            throw new ArrayIndexOutOfBoundsException("Index j cannot be less thant zero.");
        }
        if (j >= cols) {
            throw new ArrayIndexOutOfBoundsException(String.format("Index j: %d >= than number of columns: %d.", j, cols));
        }
        return data[i * cols + j];
    }

    @Override
    public Complex unsafeGet(int row, int col) {
        return data[row * cols + col];
    }

    public void set(int i, int j, Complex val) {
        // TODO clean up all these checks
        if (i < 0) {
            throw new ArrayIndexOutOfBoundsException("Index i cannot be less thant zero.");
        }
        if (i >= rows) {
            throw new ArrayIndexOutOfBoundsException(String.format("Index i: %d >= than number of rows: %d.", i, rows));
        }
        if (j < 0) {
            throw new ArrayIndexOutOfBoundsException("Index j cannot be less thant zero.");
        }
        if (j >= cols) {
            throw new ArrayIndexOutOfBoundsException(String.format("Index j: %d >= than number of columns: %d.", j, cols));
        }
        data[i * cols + j] = val;
    }

    @Override
    public void unsafeSet(int row, int col, Complex val) {
        data[row * cols + col] = val;
    }

    @Override
    public double det() {
        return 0;
    }

    @Override
    public ComplexLUDecomposition<?> LU() {
        return null;
    }

    @Override
    public ComplexQRDecomposition<?> QR() {
        return null;
    }

    public Complex[] getArray() {
        return this.data;
    }

    public Complex[] getArrayCopy() {
        return ComplexArrays.deepCopy(data);
    }
    // endregion

    // region arithmetic operations
    public ComplexMatrixDense add(ComplexMatrixDense m) {
        checkMatrixDimensions(m);
        Complex[] result = new Complex[this.rows * this.cols];
        for (int i = 0; i < this.rows * this.cols; ++i) {
            result[i] = this.data[i].add(m.data[i]);
        }
        return new ComplexMatrixDense(result, rows, cols);
    }

    public void addEquals(MatrixDense m) {
        double[] mData = m.getArray();
        checkMatrixDimensions(m);
        final int length = rows * cols;
        for (int i = 0; i < length; ++i) {
            data[i].addEquals(mData[i]);
        }
    }

    public ComplexMatrixDense subtract(MatrixDense m) {
        double[] data = m.getArray();
        checkMatrixDimensions(m);
        Complex[] result = new Complex[this.rows * this.cols];
        for (int i = 0; i < this.rows * this.cols; ++i) {
            result[i] = this.data[i].subtract(data[i]);
        }
        return new ComplexMatrixDense(result, rows, cols);
    }

    public ComplexMatrixDense subtract(ComplexMatrixDense m) {
        Complex[] data = m.getArray();
        checkMatrixDimensions(m);
        Complex[] result = new Complex[this.rows * this.cols];
        for (int i = 0; i < this.rows * this.cols; ++i) {
            result[i] = this.data[i].subtract(data[i]);
        }
        return new ComplexMatrixDense(result, rows, cols);
    }

    /**
     * Multiply a matrix by a scalar, C = s*A
     *
     * @param s scalar
     * @return s*A
     */
    public ComplexMatrixDense multiply(double s) {
        return new ComplexMatrixDense(ComplexArrays.multiplyElementWise(data, s), rows, cols);
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

    public ComplexMatrixDense multiply(MatrixDense matrix) {
        ComplexMatrixDense c = new ComplexMatrixDense(0, 0);
        multiplyOp(this, matrix, c);
        return c;
    }

    public ComplexMatrixDense multiply(ComplexMatrixDense matrix) {
        ComplexMatrixDense c = new ComplexMatrixDense(0, 0);
        multiplyOp(this, matrix, c);
        return c;
    }

    private static void multiplyOp(ComplexMatrixDense a, MatrixDense b, ComplexMatrixDense c) {
        int bRows = b.getRowCount();
        int bCols = b.getColumnCount();
        if (bRows != a.cols) {
            throw new IllegalArgumentException("Matrix inner dimensions must agree. Check that the number of" +
                    "columns of the first matrix equal the number of rows of the second matrix.");
        }
        Complex[] result = new Complex[a.rows * bCols];
        double[] bColJ = new double[a.cols];
        double[] bData = b.getArray();
        for (int j = 0; j < bCols; j++) {
            for (int k = 0; k < a.cols; k++) {
                bColJ[k] = bData[k * bCols + j];
            }
            for (int i = 0; i < a.rows; i++) {
                Complex s = new Complex();
                for (int k = 0; k < a.cols; k++) {
                    s.addEquals(a.data[k + i * a.cols].multiply(bColJ[k]));
                }
                result[i * bCols + j] = s;
            }
        }
        c.data = result;
        c.rows = a.rows;
        c.cols = bCols;
    }

    private static void multiplyOp(ComplexMatrixDense a, ComplexMatrixDense b, ComplexMatrixDense c) {
        int bRows = b.getRowCount();
        int bCols = b.getColumnCount();
        if (bRows != a.cols) {
            throw new IllegalArgumentException("Matrix inner dimensions must agree. Check that the number of" +
                    "columns of the first matrix equal the number of rows of the second matrix.");
        }
        Complex[] result = new Complex[a.rows * bCols];
        Complex[] bColJ = new Complex[a.cols];
        Complex[] bData = b.getArray();
        for (int j = 0; j < bCols; j++) {
            for (int k = 0; k < a.cols; k++) {
                bColJ[k] = bData[k * bCols + j];
            }
            for (int i = 0; i < a.rows; i++) {
                Complex s = new Complex();
                for (int k = 0; k < a.cols; k++) {
                    s.addEquals(a.data[k + i * a.cols].multiply(bColJ[k]));
                }
                result[i * bCols + j] = s;
            }
        }
        c.data = result;
        c.rows = a.rows;
        c.cols = bCols;
    }

    public void multiplyEquals(MatrixDense matrix) {
        multiplyOp(this, matrix, this);
    }

    public void multiplyEquals(ComplexMatrixDense matrix) {
        multiplyOp(this, matrix, this);
    }
    // endregion

    public boolean isEmpty() {
        if ((rows == 0 && cols == 0) || data == null || data.length == 0) {
            return true;
        }
        return false;
    }

    public ComplexMatrixDense transpose() {
        if (this.isEmpty()) {
            return new ComplexMatrixDense(null, 0, 0);
        }
        Complex[] result = new Complex[rows * cols];
        final int trows = cols;
        final int tcols = rows;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result[j * tcols + i] = data[i * cols + j];
            }
        }
        return new ComplexMatrixDense(result, trows, tcols);
    }

    public ComplexMatrixDense inv() {
        return this.solve(MatrixDense.Factory.identity(rows));
    }

    // region solve
    public ComplexMatrixDense solve(MatrixDense B) {
        return solve(fromRealMatrix(B));
    }

    public ComplexMatrixDense solve(ComplexMatrixDense B) {
        if (rows == cols) { // Matrix is Squared
            return new ComplexLUDecompositionDense(this).solve(B);
        } else if (rows > cols) { // Matrix is tall and narrow (Overdetermined system)
            return new ComplexQRDecompositionDense(this).solve(B);
        } else { // Matrix is short and wide (Under-determined system)
            throw new UnsupportedOperationException("Method not implemented yet for short and wide Matrices");
        }
    }

    // endregion

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < rows * cols; ++i) {
            if (i > 0 && i % cols == 0) {
                sb.append(System.lineSeparator());
            }
            sb.append(String.format("%s", data[i])).append(" ");
        }
        sb.setLength(sb.length() - 1);
        return sb.toString();
    }

    /**
     * Check if size(A) == size(B)
     **/
    void checkMatrixDimensions(ComplexMatrixDense B) {
        if (B.rows != rows || B.cols != cols) {
            throw new IllegalArgumentException("Matrix dimensions must agree.");
        }
    }

    /**
     * Check if size(A) == size(B)
     **/
    void checkMatrixDimensions(MatrixDense B) {
        if (B.getRowCount() != rows || B.getColumnCount() != cols) {
            throw new IllegalArgumentException("Matrix dimensions must agree.");
        }
    }

    public ComplexMatrixDense conjugateTranspose() {
        if (this.isEmpty()) {
            return new ComplexMatrixDense(null, 0, 0);
        }
        Complex[] result = new Complex[rows * cols];
        final int trows = cols;
        final int tcols = rows;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result[j * tcols + i] = data[i * cols + j].conj();
            }
        }
        return new ComplexMatrixDense(result, trows, tcols);
    }

    /***
     * Get sub-matrix
     *
     * @param row0 The initial row index.
     * @param row1 The final row index
     * @param col0 The initial column index.
     * @param col1 The final column index.
     * @return {@code A(rows(:), col0 : col1)}.
     */
    public ComplexMatrixDense subMatrix(int row0, int row1, int col0, int col1) {
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
        Complex[] data = new Complex[rowDim * colDim];
        for (int i = row0; i <= row1; ++i) {
            for (int j = col0; j <= col1; ++j) {
                data[(i - row0) * colDim + (j - col0)] = this.data[i * cols + j];
            }
        }
        return new ComplexMatrixDense(data, rowDim, colDim);
    }

    public double normFrob() {
        int i, j;
        double fac, nrm, scale;

        scale = 0.0;
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                scale = Math.max(scale,
                        Math.abs(unsafeGet(i, j).real()) + Math.abs(unsafeGet(i, j).imag()));
            }
        }
        if (scale == 0) {
            return 0.0;
        }
        if (scale < 1) {
            scale = scale * 1.0e20;
        }
        scale = 1 / scale;
        nrm = 0;
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                fac = scale * unsafeGet(i, j).real();
                nrm = nrm + fac * fac;
                fac = scale * unsafeGet(i, j).imag();
                nrm = nrm + fac * fac;
            }
        }
        return Math.sqrt(nrm) / scale;
    }

    /***
     * Get sub-matrix
     *
     * @param rows The array of row indices.
     * @param col0 The initial column index.
     * @param col1 The final column index.
     * @return {@code A(rows(:), col0 : col1)}.
     */
    public ComplexMatrixDense subMatrix(int[] rows, int col0, int col1) {
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
        Complex[] data = new Complex[rowDim * colDim];
        for (int i = 0; i < rowDim; ++i) {
            for (int j = 0; j < colDim; ++j) {
                data[i * colDim + j] = this.data[rows[i] * cols + (j + col0)];
            }
        }
        return new ComplexMatrixDense(data, rowDim, colDim);
    }

    public static final class Factory {
        private Factory() {}


        /**
         * Identity {@code Matrix}.
         *
         * @param rows The number of rows.
         * @param cols The number of columns.
         * @return {@code identity(rows, cols)}.
         */
        public static ComplexMatrixDense identity(int rows, int cols) {
            Complex[] data = new Complex[rows * cols];
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    if (i == j) {
                        data[i * cols + j] = new Complex(1, 0);
                    } else {
                        data[i * cols + j] = new Complex(0, 0);
                    }
                }
            }
            return new ComplexMatrixDense(data, rows, cols);
        }

        /**
         * Identity {@code Matrix.}
         *
         * @param n The number of rows and columns.
         * @return {@code identity(n, n)}.
         */
        public static ComplexMatrixDense identity(int n) {
            return ComplexMatrixDense.Factory.identity(n, n);
        }

        /**
         * Random {@code Matrix.}
         * @param n The number of rows and columns.
         * @return {@code random(n, n)}.
         */
        public static ComplexMatrixDense random(int n) {
            return random(n, n);
        }

        /**
         * Random {@code Matrix.}
         *
         * @param rows The number of rows
         * @param cols The number of columns.
         * @return {@code random(rows, cols)}.
         */
        public static ComplexMatrixDense random(int rows, int cols) {
            Random rand = new Random();
            Complex[] data = new Complex[rows * cols];

            for (int i = 0; i < data.length; ++i) {
                double real = rand.nextDouble() * 100.0;
                double imag = rand.nextDouble() * 100.0;
                data[i] = new Complex(real, imag);
            }
            return new ComplexMatrixDense(data, rows, cols);
        }

        public static ComplexMatrixDense zeros(int dim) {
            return new ComplexMatrixDense(dim, dim, new Complex());
        }

        public static ComplexMatrixDense zeros(int rows, int cols) {
            return new ComplexMatrixDense(rows, cols, new Complex());
        }
    }
}
