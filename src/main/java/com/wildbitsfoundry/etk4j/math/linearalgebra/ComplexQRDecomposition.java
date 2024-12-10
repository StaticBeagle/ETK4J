package com.wildbitsfoundry.etk4j.math.linearalgebra;

public abstract class ComplexQRDecomposition<T extends ComplexMatrix> {
    protected final int rows;
    protected final int cols;

    public ComplexQRDecomposition(T matrix) {
        this.rows = matrix.getRowCount();
        this.cols = matrix.getColumnCount();
    }
}
