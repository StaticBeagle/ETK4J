package com.wildbitsfoundry.etk4j.math.optimize.solvers;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.optimize.OptimizerStatusType;
import org.junit.Test;

import static org.junit.Assert.*;

public class BrentTest {

    @Test
    public void testFindRoot() {
        SolverResults<Double> b = new Brent(x ->  x * x * x - x * x + 2, -20.0, 20.0)
                .absTolerance(1e-9)
                .relTolerance(8.0 * ConstantsETK.DOUBLE_EPS)
                .iterationLimit(100)
                .solve();
        assertEquals(-0.9999999998117516, b.getValue(), 1e-12);
        assertEquals(10, b.getNumberOfIterations());
        assertEquals("Converged", b.getSolverStatus());
        assertEquals(OptimizerStatusType.CONVERGED, b.getOptimizerStatusType());
        assertEquals(5.000008185263027E-10, b.getError(), 1e-12);
        assertTrue(b.hasConverged());
    }

    @Test
    public void testFindRootExceedIterationLimit() {
        SolverResults<Double> b = new Brent(x ->  x * x * x - x * x + 2, -20.0, 20.0)
                .absTolerance(1e-9)
                .relTolerance(8.0 * ConstantsETK.DOUBLE_EPS)
                .iterationLimit(9)
                .solve();
        assertEquals(-1.0000000003117524, b.getValue(), 1e-12);
        assertEquals(10, b.getNumberOfIterations());
        assertEquals("Maximum number of iterations exceeded", b.getSolverStatus());
        assertEquals(OptimizerStatusType.MAXIMUM_NUMBER_OF_ITERATIONS_EXCEEDED, b.getOptimizerStatusType());
        assertEquals(5.000008185263027E-10, b.getError(), 1e-12);
        assertFalse(b.hasConverged());
    }
}
