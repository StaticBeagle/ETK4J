package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.util.NumArrays;
import org.junit.Test;

import java.lang.reflect.Array;
import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class QuadraticSplineZeroCenteredTest {

    @Test
    public void testNaturalSplineInterpolation() {
        double[] x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        double[] y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
        QuadraticSplineZeroCentered qs = QuadraticSplineZeroCentered.newQuadraticSpline(x, y);

        double[] expected = {0.0, 22.704, 45.408, 68.112, 90.816, 113.52000000000001, 136.224, 158.92799999999997,
                181.63199999999998, 204.33599999999998, 227.03999999999996, 250.63279999999997, 276.0032, 303.1512,
                332.0768, 362.78, 394.2363999999999, 425.4216, 456.3356, 486.97839999999997, 517.35, 549.1908000000001,
                584.2412000000002, 622.1522222222222, 660.8299999999999, 699.9255555555555, 739.438888888889,
                779.3700000000001, 819.718888888889, 860.4855555555556, 901.6700000000001};
        assertArrayEquals(expected, qs.evaluateAt(NumArrays.linSteps(0, 30)), 1e-12);

        assertEquals(394.2363999999999, qs.evaluateAt(16.0), 1e-12);
        assertEquals(31.3208, qs.differentiate(16.0), 1e-12);
        assertEquals(31.3208, qs.evaluateDerivativeAt(2, 16), 1e-12);
        System.out.println(qs);
        assertEquals(2113.5808000000075, qs.evaluateAntiDerivativeAt(2, 16), 1e-12);
        assertEquals(2969.7641333333327, qs.integrate(16), 1e-12);
        assertEquals(2969.7641333333327, qs.integrate(0.0, 16), 1e-12);
        assertEquals(1595.8758666666665, qs.integrate(11.,16.0), 1e-12);
    }

    // TODO is this below needed. If we test in the cubic spline do we have to test here as well?
//    @Test
//    public void testNaturalSplineExtrapolateLeft() {
//        CubicSpline cspline = CubicSpline.newNaturalSpline(x, y);
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_END_POINT);
//        double yi = cspline.evaluateAt(left);
//        assertEquals("Natural Spline ClampToEndPoint lower bound extrapolation", 1.3, yi, 0.0);
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_NAN);
//        yi = cspline.evaluateAt(left);
//        assertTrue("Natural Spline ClampToEndNaN lower bound extrapolation", Double.isNaN(yi));
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_ZERO);
//        yi = cspline.evaluateAt(left);
//        assertEquals("Natural Spline ClampToZero lower bound extrapolation", 0.0, yi, 0.0);
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.LINEAR);
//        yi = cspline.evaluateAt(left);
//        assertEquals("Natural Spline Linear lower bound extrapolation", 0.5474178403755874, yi, 1e-12);
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.NATURAL);
//        yi = cspline.evaluateAt(left);
//        assertEquals("Natural Spline Natural lower bound extrapolation", 1.1915492957746414, yi, 1e-12);
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.PERIODIC);
//        yi = cspline.evaluateAt(left);
//        assertEquals("Natural Spline Periodic lower bound extrapolation", 1.8500000000000005, yi, 1e-12);
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.THROW);
//        exception.expect(IndexOutOfBoundsException.class);
//        yi = cspline.evaluateAt(left);
//    }
//
//    @Test
//    public void testNaturalSplineExtrapolateRight() {
//        CubicSpline cspline = CubicSpline.newNaturalSpline(x, y);
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_END_POINT);
//        double yi = cspline.evaluateAt(right);
//        assertEquals("Natural Spline ClampToEndPoint upper bound extrapolation", 2.1, yi, 0.0);
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_NAN);
//        yi = cspline.evaluateAt(right);
//        assertTrue("Natural Spline ClampToEndNaN upper bound extrapolation", Double.isNaN(yi));
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_ZERO);
//        yi = cspline.evaluateAt(right);
//        assertEquals("Natural Spline ClampToZero upper bound extrapolation", 0.0, yi, 0.0);
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.LINEAR);
//        yi = cspline.evaluateAt(right);
//        assertEquals("Natural Spline Linear upper bound extrapolation", 3.306338028169013, yi, 1e-12);
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.NATURAL);
//        yi = cspline.evaluateAt(right);
//        assertEquals("Natural Spline Natural upper bound extrapolation", 1.6592429577464949, yi, 1e-12);
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.PERIODIC);
//        yi = cspline.evaluateAt(right);
//        assertEquals("Natural Spline Periodic upper bound extrapolation", 1.7557218309859157, yi, 1e-12);
//
//        cspline.setExtrapolationMethod(ExtrapolationMethod.THROW);
//        exception.expect(IndexOutOfBoundsException.class);
//        yi = cspline.evaluateAt(right);
//    }
}
