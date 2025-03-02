package com.wildbitsfoundry.etk4j.math.polynomial;

import com.wildbitsfoundry.etk4j.math.function.UnivariateFunction;
import org.junit.Test;

import static com.wildbitsfoundry.etk4j.math.polynomial.ChebyshevPolynomial.computeNodes;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class ChebyshevPolynomialTest {

    @Test
    public void testConstructor() {
        int n = 10; // Number of Chebyshev nodes
        double a = -1, b = 1; // Interval [a, b]
        double[] nodes = computeNodes(a, b, n);

        UnivariateFunction fn = y -> 1.0 / (1 + 25 * y * y);

        double[] y = new double[nodes.length];
        for (int i = 0; i < y.length; i++) {
            y[i] = fn.evaluateAt(nodes[i]);
        }

        double[] expected = {2.792850922403183E-17, 0.04306853365521306, -1.2244558208218224E-16, -0.09302803269526026,
                8.865098439172136E-17, 0.15787201696654918, -4.663770334977048E-19, -0.24797552395248557,
                -4.145281254381337E-17, 0.18887755738541};

        ChebyshevPolynomial poly = ChebyshevPolynomial.chebyFit(nodes, y, n - 1);
        assertArrayEquals(expected, poly.getCoefficients(), 1e-12);
    }

    @Test
    public void testEvaluateAtOneCoefficient() {
        int n = 1; // Number of Chebyshev nodes
        double a = -1, b = 1; // Interval [a, b]
        double[] nodes = computeNodes(a, b, n);

        UnivariateFunction fn = y -> 1.0 / (1 + 25 * y * y);

        double[] y = new double[nodes.length];
        for (int i = 0; i < y.length; i++) {
            y[i] = fn.evaluateAt(nodes[i]);
        }

        ChebyshevPolynomial poly = ChebyshevPolynomial.chebyFit(nodes, y, n - 1);
        assertEquals(1, poly.evaluateAt(-1), 1e-12);
        assertEquals(1, poly.evaluateAt(-0.5), 1e-12);
        assertEquals(1, poly.evaluateAt(-0.25), 1e-12);
        assertEquals(1, poly.evaluateAt(0), 1e-12);
        assertEquals(1, poly.evaluateAt(0.25), 1e-12);
        assertEquals(1, poly.evaluateAt(0.5), 1e-12);
        assertEquals(1, poly.evaluateAt(1), 1e-12);
    }

    @Test
    public void testEvaluateAtTwoCoefficient() {
        int n = 2; // Number of Chebyshev nodes
        double a = -1, b = 1; // Interval [a, b]
        double[] nodes = computeNodes(a, b, n);

        UnivariateFunction fn = y -> 1.0 / (1 + 25 * y * y);

        double[] y = new double[nodes.length];
        for (int i = 0; i < y.length; i++) {
            y[i] = fn.evaluateAt(nodes[i]);
        }

        ChebyshevPolynomial poly = ChebyshevPolynomial.chebyFit(nodes, y, n - 1);
        assertEquals(0.07407407407407411, poly.evaluateAt(-1), 1e-12);
        assertEquals(0.07407407407407411, poly.evaluateAt(-0.5), 1e-12);
        assertEquals(0.07407407407407411, poly.evaluateAt(-0.25), 1e-12);
        assertEquals(0.07407407407407411, poly.evaluateAt(0), 1e-12);
        assertEquals(0.07407407407407411, poly.evaluateAt(0.25), 1e-12);
        assertEquals(0.07407407407407411, poly.evaluateAt(0.5), 1e-12);
        assertEquals(0.07407407407407411, poly.evaluateAt(1), 1e-12);
    }

    @Test
    public void testEvaluateAt() {
        int n = 10; // Number of Chebyshev nodes
        double a = -1, b = 1; // Interval [a, b]
        double[] nodes = computeNodes(a, b, n);

        UnivariateFunction fn = y -> 1.0 / (1 + 25 * y * y);

        double[] y = new double[nodes.length];
        for (int i = 0; i < y.length; i++) {
            y[i] = fn.evaluateAt(nodes[i]);
        }

        ChebyshevPolynomial poly = ChebyshevPolynomial.chebyFit(nodes, y, n - 1);
        assertEquals(0.048814551359426445, poly.evaluateAt(-1), 1e-12);
        assertEquals(0.1193670113555115, poly.evaluateAt(-0.5), 1e-12);
        assertEquals(0.4760547551509755, poly.evaluateAt(-0.25), 1e-12);
        assertEquals(0.730821664654918, poly.evaluateAt(0), 1e-12);
        assertEquals(0.4760547551509755, poly.evaluateAt(0.25), 1e-12);
        assertEquals(0.1193670113555115, poly.evaluateAt(0.5), 1e-12);
        assertEquals(0.048814551359426445, poly.evaluateAt(1), 1e-12);
    }

    @Test
    public void testEvaluateFirstKind() {
        int n = 10; // Number of Chebyshev nodes
        double x = 0.5;
        assertEquals(-0.5, ChebyshevPolynomial.evaluateChebyshevFirstKind(n , x), 1e-12);
    }
}
