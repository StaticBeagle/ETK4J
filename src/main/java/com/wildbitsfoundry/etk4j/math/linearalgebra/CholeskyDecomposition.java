package com.wildbitsfoundry.etk4j.math.linearalgebra;

public class CholeskyDecomposition<T extends Matrix> {
    protected final int rows;
    protected final int cols;

    public CholeskyDecomposition(T matrix) {
        this.rows = matrix.getRowCount();
        this.cols = matrix.getColumnCount();
        // TODO
    }
}