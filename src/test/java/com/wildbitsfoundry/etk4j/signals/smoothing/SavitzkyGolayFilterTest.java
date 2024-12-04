package com.wildbitsfoundry.etk4j.signals.smoothing;

import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;

public class SavitzkyGolayFilterTest {

    @Test
    public void testInterpolationAtEnds() {
        double[] data = {2, 2, 5, 2, 1, 0, 1, 4, 9};
        int windowSize = 5; // Odd window size
        int polyOrder = 2; // Polynomial degree

        double[] expected = {1.6571428571428561, 3.1714285714285704, 3.542857142857143, 2.8571428571428577,
                0.6571428571428569, 0.1714285714285712, 0.9999999999999996, 3.9999999999999987, 8.999999999999998};

        double[] smoothed = new SavitzkyGolayFilter(data, windowSize, polyOrder).filter();
        assertArrayEquals(expected, smoothed, 1e-12);
    }

    @Test
    public void testMirrorAtEnds() {
        double[] data = {2, 2, 5, 2, 1, 0, 1, 4, 9};
        int windowSize = 5; // Odd window size
        int polyOrder = 2; // Polynomial degree

        double[] expected = {1.4857142857142853, 3.0285714285714285, 3.542857142857143, 2.8571428571428577,
                0.6571428571428569, 0.1714285714285712, 0.9999999999999996, 5.0285714285714285, 6.942857142857143};

        double[] smoothed = new SavitzkyGolayFilter(data, windowSize, polyOrder)
                .mode(SavitzkyGolayFilter.Mode.MIRROR)
                .filter();
        assertArrayEquals(expected, smoothed, 1e-12);
    }

    @Test
    public void testWrapAtEnds() {
        double[] data = {2, 2, 5, 2, 1, 0, 1, 4, 9};
        int windowSize = 5; // Odd window size
        int polyOrder = 2; // Polynomial degree

        double[] expected = {3.9714285714285715, 2.428571428571428, 3.542857142857143, 2.8571428571428577,
                0.6571428571428569, 0.1714285714285712, 0.9999999999999996, 5.2, 6.171428571428572};

        double[] smoothed = new SavitzkyGolayFilter(data, windowSize, polyOrder)
                .mode(SavitzkyGolayFilter.Mode.WRAP)
                .filter();
        assertArrayEquals(expected, smoothed, 1e-12);
    }

    @Test
    public void testConstantAtEnds() {
        double[] data = {2, 2, 5, 2, 1, 0, 1, 4, 9};
        int windowSize = 5; // Odd window size
        int polyOrder = 2; // Polynomial degree

        double[] expected = {1.4857142857142855, 3.1142857142857143, 3.542857142857143, 2.8571428571428577,
                0.6571428571428569, 0.1714285714285712, 0.9999999999999996, 5.2857142857142865, 5.914285714285715};

        double[] smoothed = new SavitzkyGolayFilter(data, windowSize, polyOrder)
                .mode(SavitzkyGolayFilter.Mode.CONSTANT)
                .constant(1)
                .filter();
        assertArrayEquals(expected, smoothed, 1e-12);
    }

    @Test
    public void testNearestAtEnds() {
        double[] data = {2, 2, 5, 2, 1, 0, 1, 4, 9};
        int windowSize = 5; // Odd window size
        int polyOrder = 2; // Polynomial degree

        double[] expected = {1.7428571428571429, 3.0285714285714285, 3.542857142857143, 2.8571428571428577,
                0.6571428571428569, 0.1714285714285712, 0.9999999999999996, 4.6, 7.9714285714285715};

        double[] smoothed = new SavitzkyGolayFilter(data, windowSize, polyOrder)
                .mode(SavitzkyGolayFilter.Mode.NEAREST)
                .filter();
        assertArrayEquals(expected, smoothed, 1e-12);
    }
}
