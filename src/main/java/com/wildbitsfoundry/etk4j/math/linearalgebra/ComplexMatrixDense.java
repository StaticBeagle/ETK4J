package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

import java.util.Arrays;

public class ComplexMatrixDense {
    private Complex[] data;
    private int rows;
    private int cols;

    public ComplexMatrixDense(int rows, int cols) {
        this.rows = rows;
        this.cols = cols;

        this.data = new Complex[rows * cols];
    }

    /***
     * Column packed
     * @param data
     * @param rows
     */
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

    public ComplexMatrixDense(int rows, int cols, double val) {
        this.rows = rows;
        this.cols = cols;
        data = new Complex[this.rows * this.cols];
        Arrays.fill(data, val);
    }

    /***
     * Deep copy
     * @return
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

    // region getters and setter
    public int getRowCount() {
        return rows;
    }

    public int getColumnCount() {
        return cols;
    }

    public Complex get(int i, int j) {
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

    public void set(int i, int j, Complex val) {
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

    public ComplexMatrixDense multiply(MatrixDense matrix) {
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

    public void multiplyEquals(MatrixDense matrix) {
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
        return this.solve(MatrixDense.identity(rows));
    }

    // region solve
    public ComplexMatrixDense solve(MatrixDense B) {
        return new ComplexLUDecomposition(this).solve(B);
        // Only implemented for square matrices for now.
//        if (rows == cols) { // Matrix is Squared
//            return new LUDecomposition(this).solve(B);
//        } else if (rows > cols) { // Matrix is thin (Overdetermined system)
//            return new QRDecomposition(this).solve(B);
//        } else { // Matrix is fat (Under-determined system)
//            QRDecomposition qr = this.transpose().QR();
//            Matrix R1 = fwdSubsSolve(qr.getRT(), B);
//            R1.appendRows(cols - R1.rows);
//            return qr.QmultiplyX(R1);
//        }
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

    public static void main(String[] args) {
        MatrixDense A = new MatrixDense(new double[][]{{-2, -1}, {1, 0}});
        MatrixDense B = new MatrixDense(new double[][]{{1}, {0}});
        MatrixDense C = new MatrixDense(new double[][]{{1, 2}});
        MatrixDense D = new MatrixDense(new double[][]{{1}});

        double w = 100;
        Complex jw = Complex.fromImaginary(w);
        ComplexMatrixDense gg = MatrixDense.identity(A.getRowCount()).multiply(jw).subtract(A);
        ComplexMatrixDense inv = gg.inv();
        System.out.println(gg);
        System.out.println(inv);
    }
}
