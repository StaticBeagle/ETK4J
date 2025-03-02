package com.wildbitsfoundry.etk4j.math.optimize.minimizer;

import com.wildbitsfoundry.etk4j.math.optimize.OptimizerStatusType;
import org.junit.Test;

import java.util.function.BiFunction;

import static org.junit.Assert.*;

public class BrentTest {

    @Test
    public void testFindMinimum() {
        BiFunction<Double, Object[], Double> function = (x, obj) -> Math.pow(x, 4) + 20 * Math.pow(x, 3) - x * x + 2;
        MinimizerResults<Double> b = new Brent(function, -20, 20)
                .tolerance(1.48e-8)
                .iterationLimit(500)
                .minimize();
        assertEquals(-15.033259534575702, b.getValue(), 1e-12);
        assertEquals(-17098.49963085873, b.getFunctionValue(), 1e-12);
        assertEquals(13, b.getNumberOfIterations());
        assertEquals("Converged", b.getMinimizerStatus());
        assertEquals(OptimizerStatusType.CONVERGED, b.getOptimizerStatusType());
        assertTrue(b.hasConverged());
    }

    @Test
    public void testFindMinimumExceedIterationLimit() {
        BiFunction<Double, Object[], Double> function = (x, obj) -> Math.pow(x, 4) + 20 * Math.pow(x, 3) - x * x + 2;
        MinimizerResults<Double> b = new Brent(function, -20, 20)
                .tolerance(1.48e-8)
                .iterationLimit(10)
                .minimize();
        assertEquals(-15.033230701070634, b.getValue(), 1e-12);
        assertEquals(-17098.49963048077, b.getFunctionValue(), 1e-12);
        assertEquals(10, b.getNumberOfIterations());
        assertEquals("Maximum number of iterations exceeded", b.getMinimizerStatus());
        assertEquals(OptimizerStatusType.MAXIMUM_NUMBER_OF_ITERATIONS_EXCEEDED, b.getOptimizerStatusType());
        assertFalse(b.hasConverged());
    }
}
