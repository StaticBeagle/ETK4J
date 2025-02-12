package com.wildbitsfoundry.etk4j.math.linearalgebra;
// TODO do we need more (abstract) methods here?
public abstract class QRDecomposition<T extends Matrix> {
    protected final int rows;
    protected final int cols;

    public QRDecomposition(T matrix) {
        this.rows = matrix.getRowCount();
        this.cols = matrix.getColumnCount();
    }
}
