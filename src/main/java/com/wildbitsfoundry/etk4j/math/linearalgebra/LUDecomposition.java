package com.wildbitsfoundry.etk4j.math.linearalgebra;

public abstract class LUDecomposition <T extends Matrix> {

    protected final int rows;
    protected final int cols;

    public LUDecomposition(T matrix) {
        this.rows = matrix.getRowCount();
        this.cols = matrix.getColumnCount();
    }

    public abstract boolean isSingular();

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
