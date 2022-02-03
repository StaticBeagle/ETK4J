package com.wildbitsfoundry.etk4j.math.optimize.solvers;

import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class NewtonRaphsonTest {

    @Test
    public void testVanillaNewtonRaphsonMethod() {
        SolverResults<Double> nr = new NewtonRaphson(x ->  x * x * x - x * x + 2, -20.0)
                .derivative(x -> 3 * x * x - 2 * x)
                .absTolerance(1e-9)
                .relTolerance(0.0)
                .iterationLimit(100)
                .solve();
        assertEquals(-1.0, nr.getValue(), 1e-12);
        assertEquals(12, nr.getNumberOfIterations());
        assertEquals("Converged", nr.getSolverStatus());
        assertEquals(1.1102230246251565E-14, nr.getError(), 1e-12);
        assertTrue(nr.hasConverged());
    }

    @Test
    public void testHalleysMethod() {
        SolverResults<Double> nr = new NewtonRaphson(x ->  x * x * x - x * x + 2, -20.0)
                .derivative(x -> 3 * x * x - 2 * x)
                .secondDerivative(x -> 6 * x)
                .solve();
        assertEquals(-1.0, nr.getValue(), 1e-12);
        assertEquals(9, nr.getNumberOfIterations());
        assertEquals("Converged", nr.getSolverStatus());
        assertEquals(0.0, nr.getError(), 1e-12);
        assertTrue(nr.hasConverged());
    }

    @Test
    public void testSecantMethod() {
        SolverResults<Double> nr = new NewtonRaphson(x ->  x * x * x - x * x + 2, -20.0)
                .solve();
        assertEquals(-1.0000000000000002, nr.getValue(), 1e-12);
        assertEquals(16, nr.getNumberOfIterations());
        assertEquals("Converged", nr.getSolverStatus());
        assertEquals(2.2784485453897219E-7, nr.getError(), 1e-12);
        assertTrue(nr.hasConverged());
    }

    @Test
    public void testSecantMethodWithSecondGess() {
        SolverResults<Double> nr = new NewtonRaphson(x ->  x * x * x - x * x + 2, -20.0)
                .secondInitialGuess(-19.0)
                .solve();
        assertEquals(-1.0, nr.getValue(), 1e-12);
        assertEquals(16, nr.getNumberOfIterations());
        assertEquals("Converged", nr.getSolverStatus());
        assertEquals(1.335008938951887E-7, nr.getError(), 1e-12);
        assertTrue(nr.hasConverged());
    }
}
