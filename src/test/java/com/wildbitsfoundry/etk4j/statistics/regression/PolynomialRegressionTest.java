package com.wildbitsfoundry.etk4j.statistics.regression;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class PolynomialRegressionTest {
    PolynomialRegression pr;

    @Before
    public void before() {
        double[] x = { 1, 2, 3, 4 };
        double[] y = { 1, 3.5, 8.4, 17.2 };
        pr = new PolynomialRegression(x, y, 2);
    }

    @Test
    public void testBeta() {
        double[] beta = { 1.575, -2.525000000000001, 2.025000000000002 };
        assertArrayEquals(beta, pr.beta(), 1e-12);
    }

    @Test
    public void testResiduals() {
        double[] residuals = { -0.07500000000000129, 0.22499999999999964, -0.22499999999999787, 0.07499999999999929 };
        assertArrayEquals(residuals, pr.residuals(), 1e-12);
    }

    @Test
    public void testNormOfResiduals() {
        assertEquals(0.33541019662496696, pr.normOfResiduals(), 1e-12);
    }

    @Test
    public void testR() {
        assertEquals(0.999632639553955, pr.R(), 1e-12);
    }

    @Test
    public void testR2() {
        assertEquals(0.9992654140616073, pr.R2(), 1e-12);
    }

    @Test
    public void testSSE() {
        assertEquals(0.112499999999999, pr.SSE(), 1e-12);
    }

    @Test
    public void testSSR() {
        assertEquals(153.14749999999998, pr.SSR(), 1e-12);
    }
}
