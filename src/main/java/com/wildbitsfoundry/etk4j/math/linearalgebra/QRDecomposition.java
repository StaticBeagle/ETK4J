package com.wildbitsfoundry.etk4j.math.linearalgebra;

public abstract class QRDecomposition<T extends Matrix> {
    protected final int rows;
    protected final int cols;

    public QRDecomposition(T matrix) {
        this.rows = matrix.getRowCount();
        this.cols = matrix.getColumnCount();
    }
}
