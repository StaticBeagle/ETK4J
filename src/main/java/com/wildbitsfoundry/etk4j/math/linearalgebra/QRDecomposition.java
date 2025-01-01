package com.wildbitsfoundry.etk4j.math.linearalgebra;
// TODO do we need more (abstract) methods here?
public abstract class QRDecomposition<T extends Matrix> {
    protected final int rows;
    protected final int cols;
    // TODO decide on what to call the V vectors matrix and also compact Q and compact R getter interface
    public QRDecomposition(T matrix) {
        this.rows = matrix.getRowCount();
        this.cols = matrix.getColumnCount();
    }
}
