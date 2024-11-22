package com.wildbitsfoundry.etk4j.math.linearalgebra;

import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;

public class JacobiMethodSolverTest {

    @Test
    public void testSolve() {
        double[][] A = {{4, -1, 0, 1, 0},
                {-1, 4, -1, 0, 1},
                {0, -1, 4, -1, 0},
                {1, 0, -1, 4, -1},
                {0, 1, 0, -1, 4}};

        double[] b = {100, 100, 100, 100, 100};

        double[] x = new JacobiMethodSolver(new MatrixDense(A), b)
                .iterationLimit(20)
                .tolerance(0.000001)
                .initialGuess(new double[] {0, 0, 0, 0, 0})
                .solve().getValue();
        assertArrayEquals(new double[] {25.0, 35.71428544819355, 42.85714253783226, 35.71428544819355, 25.0}, x, 1e-12);
    }
}
