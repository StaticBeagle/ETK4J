package com.wildbitsfoundry.etk4j.math.polynomial;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class LagrangePolynomialTest {

    @Test
    public void testEvaluateAt() {
        double[] x = {1, 2, 3, 4};
        double[] y = {1, 10, 12, 15};
        LagrangePolynomial lagrangePolynomial = new LagrangePolynomial(x, y);
        assertEquals(1, lagrangePolynomial.evaluateAt(1), 1e-12);
        assertEquals(6.875000000000001, lagrangePolynomial.evaluateAt(1.5), 1e-12);
        assertEquals(10, lagrangePolynomial.evaluateAt(2), 1e-12);
        assertEquals(11.374999999999993, lagrangePolynomial.evaluateAt(2.5), 1e-12);
        assertEquals(12, lagrangePolynomial.evaluateAt(3), 1e-12);
        assertEquals(12.875000000000002, lagrangePolynomial.evaluateAt(3.5), 1e-12);
        assertEquals(15, lagrangePolynomial.evaluateAt(4), 1e-12);
    }

    @Test
    public void testPolyFit() {
        double[] x = {1, 2, 3, 4};
        double[] y = {1, 10, 12, 15};
        LagrangePolynomial lagrangePolynomial = LagrangePolynomial.lagrangeFit(x, y);
        assertEquals(1, lagrangePolynomial.evaluateAt(1), 1e-12);
        assertEquals(6.875000000000001, lagrangePolynomial.evaluateAt(1.5), 1e-12);
        assertEquals(10, lagrangePolynomial.evaluateAt(2), 1e-12);
        assertEquals(11.374999999999993, lagrangePolynomial.evaluateAt(2.5), 1e-12);
        assertEquals(12, lagrangePolynomial.evaluateAt(3), 1e-12);
        assertEquals(12.875000000000002, lagrangePolynomial.evaluateAt(3.5), 1e-12);
        assertEquals(15, lagrangePolynomial.evaluateAt(4), 1e-12);
    }
}
