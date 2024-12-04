package com.wildbitsfoundry.etk4j.math.linearalgebra;

import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;

public class GaussSeidelSolverTest {

    @Test
    public void testSolve() {
        double[][] A = {{4, -1, 0, 1, 0},
                {-1, 4, -1, 0, 1},
                {0, -1, 4, -1, 0},
                {1, 0, -1, 4, -1},
                {0, 1, 0, -1, 4}};

        double[] b = {100, 100, 100, 100, 100};

        double[] x = new GaussSeidelSolver(new MatrixDense(A), b)
                .iterationLimit(20)
                .tolerance(0.000001)
                .initialGuess(new double[] {0, 0, 0, 0, 0})
                .solve().getValue();
        assertArrayEquals(new double[] {25.000000018889274, 35.71428570361077, 42.857142835082136, 35.7142856993259, 24.999999998928782}, x, 1e-12);
    }
}
