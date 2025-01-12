package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public abstract class ComplexMatrix {
    protected int rows;
    protected int cols;

    public ComplexMatrix() {}

    /**
     * Constructs a {@code Matrix}.
     *
     * @param rows The number of rows of the {@code Matrix}..
     * @param cols The number of columns of the {@code Matrix}..
     */
    public ComplexMatrix(int rows, int cols) {
        this.rows = rows;
        this.cols = cols;
    }

    /***
     * Deep copy.
     * @return A deep copy of the {@code Matrix}.
     */
    public abstract ComplexMatrix copy();

    /**
     * Retrieve value from {@code Matrix} at a given position.
     *
     * @param row The row index.
     * @param col The column index.
     * @return The value at {@code A(i, j)}.
     */
    public abstract Complex get(int row, int col);

    /**
     * Retrieve value from {@code Matrix} at a given position without any boundary checks.
     *
     * @param row The row index.
     * @param col The column index.
     * @return The value at {@code A(i, j)}.
     */
    public abstract Complex unsafeGet(int row, int col);

    /**
     * Set the value of the {@code Matrix} at a given position.
     *
     * @param row   The row index.
     * @param col   The column index.
     * @param val The value used to set {@code A(i, j) = val}.
     */
    public abstract void set(int row, int col, Complex val);

    /**
     * Set the value of the {@code Matrix} at a given position without any boundary checks.
     *
     * @param row   The row index.
     * @param col   The column index.
     * @param val The value used to set {@code A(i, j) = val}.
     */
    public abstract void unsafeSet(int row, int col, Complex val);

    /**
     * {@code Matrix} determinant.
     *
     * @return The determinant of the squared {@code Matrix}.
     */
    public abstract double det();

//    /**
//     * Matrix trace.
//     *
//     * @return sum of the diagonal elements.
//     */
//    public double trace() {
//        double t = 0;
//        for (int i = 0; i < Math.min(rows, cols); ++i) {
//            t += unsafeGet(i, i);
//        }
//        return t;
//    }

    /**
     * Is the {@code Matrix} squared.
     *
     * @return {@code true} if the matrix is squared (same number of rows and columns) or {@code false} otherwise.
     */
    public boolean isSquared() {
        return rows == cols;
    }

    // region decompositions

    /**
     * LU decomposition of the {@code Matrix}.
     *
     * @return The {@link LUDecompositionDense} of the {@code Matrix}.
     * @see <a href="https://en.wikipedia.org/wiki/LU_decomposition">LU Decomposition</a>
     */
    public abstract ComplexLUDecomposition<?> LU();

    /**
     * QR decomposition of the {@code Matrix}.
     *
     * @return The {@link QRDecompositionDense} of the {@code Matrix}.
     * @see <a href="https://en.wikipedia.org/wiki/QR_decomposition">QR Decomposition</a>
     */
    public abstract ComplexQRDecomposition<?> QR();
    /**
     * Cholesky decomposition of the {@code Matrix}.
     *
     * @return The {@link CholeskyDecompositionDense} of the {@code Matrix}.
     * @see <a href="https://en.wikipedia.org/wiki/Cholesky_decomposition">Cholesky Decomposition</a>
     */
    //public abstract ComplexCholeskyDecomposition<?> Chol();
    // endregion

    /***
     * Get matrix diagonal
     *
     * @return Array containing the diagonal of the matrix
     */
    public Complex[] diag() {
        if (!this.isSquared()) {
            throw new NonSquareMatrixException("Matrix is not squared");
        }
        final int dim = rows;
        Complex[] diag = new Complex[dim];
        for (int i = 0; i < dim; ++i) {
            diag[i] = this.unsafeGet(i, i);
        }
        return diag;
    }

    /**
     * Is the {@code Matrix} empty.
     *
     * @return {@code true} if the number of columns and rows are equal to zero or the internal data is null or the
     * length of the internal data is equal to zero. Returns {@code false} otherwise.
     */
    public abstract boolean isEmpty();

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

}
