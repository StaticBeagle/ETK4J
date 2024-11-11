package com.wildbitsfoundry.etk4j.math.linearalgebra;

import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;

public class MatricesTest {

    @Test
    public void testForwardSubstitutionSolve() {
        MatrixDense L = MatrixDense.magic(3).LU().getL();
        MatrixDense b = new MatrixDense(new double[][] {{1}, {2}, {3}});
        assertArrayEquals(L.solve(b).getArray(), Matrices.forwardSubstitutionSolve(L, b).getArray(), 1e-12);
    }
}
