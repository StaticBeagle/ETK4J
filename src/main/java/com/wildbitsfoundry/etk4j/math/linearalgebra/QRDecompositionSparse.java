package com.wildbitsfoundry.etk4j.math.linearalgebra;

public class QRDecompositionSparse extends QRDecomposition<MatrixSparse> {
	protected double[] _data;
	protected double[] _rdiag;

	public QRDecompositionSparse(MatrixSparse matrix) {
		super(matrix);
		// TODO
	}
}
