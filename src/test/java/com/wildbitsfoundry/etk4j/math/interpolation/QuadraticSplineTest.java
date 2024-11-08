package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.util.DoubleArrays;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.*;

public class QuadraticSplineTest {

    @Test
    public void testNaturalSplineInterpolation() {
        double[] x = {0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        double[] y = {0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
        QuadraticSpline qs = QuadraticSpline.newNaturalSpline(x, y);

        double[] expected = {0.0, 24.270666666666678, 48.1931851851852, 71.76755555555557, 94.99377777777781,
                117.87185185185187, 140.4017777777778, 162.5835555555556, 184.4171851851852, 205.90266666666668,
                227.04, 249.24020740740738, 273.9143111111111, 301.06231111111106, 330.68420740740737, 362.78,
                395.62899259259257, 427.51048888888886, 458.42448888888885, 488.3709925925926, 517.35,
                548.1463555555556, 583.5449037037038, 622.8833333333333, 662.71, 702.5366666666666,
                742.3633333333333, 782.19, 822.0166666666667, 861.8433333333332, 901.67};
        assertArrayEquals(expected, qs.evaluateAt(DoubleArrays.linSteps(0, 30)), 1e-12);

        assertEquals(395.62899259259257, qs.evaluateAt(16.0), 1e-12);
        assertEquals(32.36524444444445, qs.differentiate(16.0), 1e-12);
        assertEquals(32.36524444444445, qs.evaluateDerivativeAt(2, 16), 1e-12);
        assertEquals(379.2851209876543, qs.evaluateAntiDerivativeAt(2, 16), 1e-12);
        assertEquals(2992.277713580247, qs.integrate(16), 1e-12);
        assertEquals(2992.277713580247, qs.integrate(0.0, 16), 1e-12);
        assertEquals(1590.1314222222222, qs.integrate(11., 16.0), 1e-12);
    }

    @Test
    public void testClampedSplineInterpolation() {
        double[] x = {0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        double[] y = {0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
        QuadraticSpline qs = QuadraticSpline.newClampedSpline(x, y, 1, 3);

        System.out.println(Arrays.toString(qs.evaluateAt(DoubleArrays.linSteps(0, 30))));

        double[] expected = {0.0, 30.214199999999995, 58.75946666666666, 85.63579999999999, 110.84319999999998,
                134.38166666666663, 156.25119999999998, 176.45179999999996, 194.98346666666663, 211.84619999999998,
                227.04, 243.95706666666666, 265.9896, 293.1376, 325.4010666666667, 362.78, 400.9121333333333,
                435.43519999999995, 466.3492, 493.65413333333333, 517.35, 544.1840000000001, 580.9033333333334,
                640.0691111111112, 706.902, 763.9144444444444, 811.1064444444444, 848.4780000000001, 876.0291111111112,
                893.7597777777778, 901.6700000000001};
        assertArrayEquals(expected, qs.evaluateAt(DoubleArrays.linSteps(0, 30)), 1e-12);

        assertEquals(400.9121333333333, qs.evaluateAt(16.0), 1e-12);
        assertEquals(36.3276, qs.differentiate(16.0), 1e-12);
        assertEquals(36.3276, qs.evaluateDerivativeAt(2, 16), 1e-12);
        assertEquals(382.14682222222217, qs.evaluateAntiDerivativeAt(2, 16), 1e-12);
        assertEquals(3077.688488888888, qs.integrate(16), 1e-12);
        assertEquals(3077.688488888888, qs.integrate(0.0, 16), 1e-12);
        assertEquals(1568.3384666666661, qs.integrate(11., 16.0), 1e-12);
    }
}
