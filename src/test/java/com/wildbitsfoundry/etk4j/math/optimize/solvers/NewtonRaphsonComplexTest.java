package com.wildbitsfoundry.etk4j.math.optimize.solvers;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.optimize.OptimizerStatusType;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class NewtonRaphsonComplexTest {

    @Test
    public void testVanillaNewtonRaphsonMethod() {
        SolverResults<Complex> nr = new NewtonRaphsonComplex(z -> z.pow(2).add(1.0), new Complex(5, -2))
                .derivative(z -> z.multiply(2.0))
                .absTolerance(1e-9)
                .relTolerance(0.0)
                .iterationLimit(100)
                .solve();
        assertEquals(Complex.fromReal(0.0).real(), nr.getValue().real(), 1e-12);
        assertEquals(Complex.fromImaginary(-1.0).imag(), nr.getValue().imag(), 1e-12);
        assertEquals(9, nr.getNumberOfIterations());
        assertEquals("Converged", nr.getSolverStatus());
        assertEquals(OptimizerStatusType.CONVERGED, nr.getOptimizerStatusType());
        assertEquals(2.472314700550465E-15, nr.getError(), 1e-12);
        assertTrue(nr.hasConverged());
    }

    @Test
    public void testHalleysMethod() {
        SolverResults<Complex> nr = new NewtonRaphsonComplex(z -> z.pow(2).add(1.0), new Complex(5, -2))
                .derivative(z -> z.multiply(2.0))
                .secondDerivative(z -> Complex.fromReal(2.0))
                .solve();
        assertEquals(Complex.fromReal(0.0).real(), nr.getValue().real(), 1e-12);
        assertEquals(Complex.fromImaginary(-1.0).imag(), nr.getValue().imag(), 1e-12);
        assertEquals(6, nr.getNumberOfIterations());
        assertEquals("Converged", nr.getSolverStatus());
        assertEquals(OptimizerStatusType.CONVERGED, nr.getOptimizerStatusType());
        assertEquals(1.3892028701619707E-14, nr.getError(), 1e-12);
        assertTrue(nr.hasConverged());
    }

    @Test
    public void testSecantMethod() {
        SolverResults<Complex> nr = new NewtonRaphsonComplex(z -> z.pow(2).add(1.0), new Complex(5, -2))
                .solve();
        assertEquals(Complex.fromReal(0.0).real(), nr.getValue().real(), 1e-12);
        assertEquals(Complex.fromImaginary(-1.0).imag(), nr.getValue().imag(), 1e-12);
        assertEquals(12, nr.getNumberOfIterations());
        assertEquals("Converged", nr.getSolverStatus());
        assertEquals(OptimizerStatusType.CONVERGED, nr.getOptimizerStatusType());
        assertEquals(8.19014815240895E-9, nr.getError(), 1e-12);
        assertTrue(nr.hasConverged());
    }

    @Test
    public void testSecantMethodWithSecondInitialGuess() {
        SolverResults<Complex> nr = new NewtonRaphsonComplex(z -> z.pow(2).add(1.0), new Complex(5, -2))
                .secondInitialGuess(new Complex(4, -1.5))
                .solve();
        assertEquals(Complex.fromReal(0.0).real(), nr.getValue().real(), 1e-12);
        assertEquals(Complex.fromImaginary(-1.0).imag(), nr.getValue().imag(), 1e-12);
        assertEquals(12, nr.getNumberOfIterations());
        assertEquals("Converged", nr.getSolverStatus());
        assertEquals(OptimizerStatusType.CONVERGED, nr.getOptimizerStatusType());
        assertEquals(2.309006245526933E-9, nr.getError(), 1e-12);
        assertTrue(nr.hasConverged());
    }
}
