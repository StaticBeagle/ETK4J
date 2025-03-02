package com.wildbitsfoundry.etk4j.math.optimize.solver;

import com.wildbitsfoundry.etk4j.constant.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.optimize.OptimizerStatusType;
import org.junit.Test;

import static org.junit.Assert.*;

public class RidderTest {

    @Test
    public void testFindRoot() {
        SolverResults<Double> r = new Ridder(x ->  x * x * x - x * x + 2, -20.0, 20.0)
                .absTolerance(1e-9)
                .relTolerance(8.0 * ConstantsETK.DOUBLE_EPS)
                .iterationLimit(10)
                .solve();
        assertEquals(-1.0, r.getValue(), 1e-12);
        assertEquals(10, r.getNumberOfIterations());
        assertEquals("Converged", r.getSolverStatus());
        assertEquals(OptimizerStatusType.CONVERGED, r.getOptimizerStatusType());
        assertEquals(2.968403300940281E-12, r.getError(), 1e-12);
        assertTrue(r.hasConverged());
    }

    @Test
    public void testFindRootExceedIterationLimit() {
        SolverResults<Double> r = new Ridder(x ->  x * x * x - x * x + 2, -20.0, 20.0)
                .absTolerance(1e-9)
                .relTolerance(8.0 * ConstantsETK.DOUBLE_EPS)
                .iterationLimit(9)
                .solve();
        assertEquals(-0.999999999997031, r.getValue(), 1e-12);
        assertEquals(10, r.getNumberOfIterations());
        assertEquals("Maximum number of iterations exceeded", r.getSolverStatus());
        assertEquals(OptimizerStatusType.MAXIMUM_NUMBER_OF_ITERATIONS_EXCEEDED, r.getOptimizerStatusType());
        assertEquals(7.21692272609431E-9, r.getError(), 1e-12);
        assertFalse(r.hasConverged());
    }
}
