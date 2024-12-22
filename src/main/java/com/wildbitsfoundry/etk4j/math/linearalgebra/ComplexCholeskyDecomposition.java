package com.wildbitsfoundry.etk4j.math.linearalgebra;

public class ComplexCholeskyDecomposition<T extends ComplexMatrix> {
    protected final int rows;
    protected final int cols;

    public ComplexCholeskyDecomposition(T matrix) {
        this.rows = matrix.getRowCount();
        this.cols = matrix.getColumnCount();
    }
}
