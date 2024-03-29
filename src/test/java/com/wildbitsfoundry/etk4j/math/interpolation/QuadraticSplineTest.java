package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.util.DoubleArrays;
import org.junit.Test;

import static org.junit.Assert.*;

public class QuadraticSplineTest {

    @Test
    public void testNaturalSplineInterpolation() {
        double[] x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        double[] y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
        QuadraticSpline qs = QuadraticSpline.newNaturalSpline(x, y);

        double[] expected = {0.0, 2.2704, 9.0816, 20.4336, 36.3264, 56.760000000000005, 81.7344, 111.2496, 145.3056,
                183.9024, 227.04, 268.796, 303.248, 330.39599999999996, 350.24, 362.78, 376.0732, 398.17679999999996,
                429.09079999999994, 468.8152, 517.35, 562.8132, 593.3228, 611.5570222222223, 633.5852,
                662.0855555555555, 697.0580888888888, 738.5028, 786.4196888888889, 840.8087555555555, 901.67};
        assertArrayEquals(expected, qs.evaluateAt(DoubleArrays.linSteps(0, 30)), 1e-12);

        assertEquals(376.0732, qs.evaluateAt(16.0), 1e-12);
        assertEquals(17.6984, qs.differentiate(16.0), 1e-12);
        assertEquals(17.6984, qs.evaluateDerivativeAt(2, 16), 1e-12);
        assertEquals(368.69239999999996, qs.evaluateAntiDerivativeAt(2, 16), 1e-12);
        assertEquals(2676.1257333333333, qs.integrate(16), 1e-12);
        assertEquals(2676.1257333333333, qs.integrate(0.0, 16), 1e-12);
        assertEquals(1670.7990666666667, qs.integrate(11.,16.0), 1e-12);
    }

    @Test
    public void testClampedSplineInterpolation() {
        double[] x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        double[] y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
        QuadraticSpline qs = QuadraticSpline.newClampedSpline(x, y, 1);

        double[] expected = {0.0, 13.704, 29.408, 47.112, 66.816, 88.52000000000001, 112.224, 137.928, 165.632, 195.336,
                227.04, 258.6328, 288.0032, 315.1512, 340.0768, 362.78, 386.23639999999995, 413.42159999999996,
                444.3356, 478.97839999999997, 517.35, 555.1908000000001, 588.2412, 617.4855555555556, 648.83,
                683.2588888888889, 720.7722222222222, 761.37, 805.0522222222222, 851.8188888888889, 901.67};
        assertArrayEquals(expected, qs.evaluateAt(DoubleArrays.linSteps(0, 30)), 1e-12);

        assertEquals(386.23639999999995, qs.evaluateAt(16.0), 1e-12);
        assertEquals(25.3208, qs.differentiate(16.0), 1e-12);
        assertEquals(25.3208, qs.evaluateDerivativeAt(2, 16), 1e-12);
        assertEquals(374.1974666666666, qs.evaluateAntiDerivativeAt(2, 16), 1e-12);
        assertEquals(2840.4307999999996, qs.integrate(16), 1e-12);
        assertEquals(2840.4307999999996, qs.integrate(0.0, 16), 1e-12);
        assertEquals(1628.8758666666663, qs.integrate(11.,16.0), 1e-12);
    }
}
