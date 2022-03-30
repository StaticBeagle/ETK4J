package com.wildbitsfoundry.etk4j.math.optimize.minimizers;

import org.junit.Test;

import java.util.function.BiFunction;

import static org.junit.Assert.*;

public class GoldenSectionTest {

    @Test
    public void testFindMinimum() {
        BiFunction<Double, Object[], Double> function = (x, obj) -> Math.pow(x, 4) + 20 * Math.pow(x, 3) - x * x + 2;
        MinimizerResults<Double> b = new GoldenSection(function, -20, 20)
                .tolerance(1.48e-8)
                .iterationLimit(500)
                .minimize();
        assertEquals(-15.033259460970806, b.getValue(), 1e-12);
        assertEquals(-17098.499630858736, b.getFunctionValue(), 1e-12);
        assertEquals(47, b.getNumberOfIterations());
        assertEquals("Converged", b.getMinimizerStatus());
        assertTrue(b.hasConverged());
    }

    @Test
    public void testFindMinimumExceedIterationLimit() {
        BiFunction<Double, Object[], Double> function = (x, obj) -> Math.pow(x, 4) + 20 * Math.pow(x, 3) - x * x + 2;
        MinimizerResults<Double> b = new GoldenSection(function, -20, 20)
                .tolerance(1.48e-8)
                .iterationLimit(10)
                .minimize();
        assertEquals(-15.01552810007571, b.getValue(), 1e-12);
        assertEquals(-17098.35742944805, b.getFunctionValue(), 1e-12);
        assertEquals(10, b.getNumberOfIterations());
        assertEquals("Maximum number of iterations exceeded", b.getMinimizerStatus());
        assertFalse(b.hasConverged());
    }
}
