package com.wildbitsfoundry.etk4j.math.specialfunctions;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class EllipticTest {

    @Test
    public void testCompleteOfFirstKind() {
        double k = 0.1;
        double expected = 1.6124413487202198;
        double y = Elliptic.completeEllipticIntegralFirstKind(k);

        assertEquals(expected, y, 1e-12);
    }

    @Test
    public void testIncompleteOfFirstKind() {
        double expected = 1.2261911708835167;
        double y = Elliptic.incompleteEllipticIntegralFirstKind(1.0, 1.0);

        assertEquals(expected, y, 1e-12);
    }

    @Test
    public void testCompleteOfSecondKind() {
        double k = 0.1;
        double expected = 1.5307576368977631;
        double y = Elliptic.completeEllipticIntegralSecondKind(k);

        assertEquals(expected, y, 1e-12);
    }

    @Test
    public void testIncompleteOfSecondKind() {
        double expected = 0.8414709848078963;
        double y = Elliptic.incompleteEllipticIntegralSecondKind(1.0, 1.0);

        assertEquals(expected, y, 1e-12);
    }
}
