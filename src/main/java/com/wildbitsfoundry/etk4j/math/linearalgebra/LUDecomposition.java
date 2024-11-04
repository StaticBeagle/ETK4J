package com.wildbitsfoundry.etk4j.math.linearalgebra;

public abstract class LUDecomposition <T extends Matrix> {

    protected final int rows;
    protected final int cols;

    public LUDecomposition(T matrix) {
        this.rows = matrix.getRowCount();
        this.cols = matrix.getColumnCount();
    }

    public abstract boolean isNonSingular();

    /**
     * Return lower triangular factor
     *
     * @return L
     */

    public abstract T getL();

    /**
     * Return upper triangular factor
     *
     * @return U
     */

    public abstract T getU();

    /**
     * Return pivot permutation vector
     *
     * @return piv
     */

//    public int[] getPivot() {
//        return Arrays.copyOf(_pivot, _pivot.length);
//    }
//
//    /**
//     * Return pivot permutation vector as a one-dimensional double array
//     *
//     * @return (double) piv
//     */
//
//    public double[] getPivotAsDouble() {
//        final int rows = _rows;
//        double[] vals = new double[rows];
//        for (int i = 0; i < rows; i++) {
//            vals[i] = (double) _pivot[i];
//        }
//        return vals;
//    }
//
//    public double det() {
//        if (_rows != _cols) {
//            throw new IllegalArgumentException("Matrix must be square.");
//        }
//        double det = _pivotsign;
//        for (int i = 0; i < _cols; i++) {
//            det *= this._data[i * _cols + i];
//        }
//        return det;
//    }

    /**
     * Solve A*X = B
     *
     * @param b
     *            A Matrix with as many rows as A and any number of columns.
     * @return X so that L*U*X = B(piv,:)
     * @exception IllegalArgumentException
     *                Matrix row dimensions must agree.
     * @exception RuntimeException
     *                Matrix is singular.
     */

    public abstract T solve(double[] b);
}
