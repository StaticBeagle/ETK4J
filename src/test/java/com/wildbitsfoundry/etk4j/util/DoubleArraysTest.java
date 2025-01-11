package com.wildbitsfoundry.etk4j.util;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class DoubleArraysTest {


    @Test
    public void testLinspace() {
        double[] a = new double[]{-1.0, -0.5, 0.0, 0.5, 1.0};
        assertArrayEquals(a, DoubleArrays.linSpace(-1.0, 1.0, 5), 1e-12);
    }

    @Test
    public void testLinspaceStep() {
        double[] a = new double[]{-1.0, -0.5, 0.0, 0.5, 1.0};
        assertArrayEquals(a, DoubleArrays.linSteps(-1.0, 1.0, 0.5), 1e-12);
    }

    @Test
    public void testLogspace() {
        double[] b = new double[]{0.1000, 0.31622776601683794, 1.0000, 3.1622776601683795, 10.0000};
        assertArrayEquals(b, DoubleArrays.logSpace(-1, 1, 5), 1e-12);
    }

    @Test
    public void testConv() {
        double[] a = new double[]{1, 1, 2, 3, 5, 8, 13, 21};
        double[] b = new double[]{5, 6};
        double[] c = new double[]{0, 1, 1, 0, 2, 3, 5, 0, 8, 13, 21};

        double[] conv = {1.0, 2.0, 5.0, 10.0, 20.0, 38.0, 71.0, 130.0, 167.0, 242.0, 320.0, 418.0, 505.0, 546.0, 441.0};
        assertArrayEquals(conv, DoubleArrays.convolve(a, a), 1e-12);

        conv = new double[]{0.0, 1.0, 2.0, 3.0, 7.0, 13.0, 25.0, 38.0, 71.0, 88.0, 125.0, 192.0,
                249.0, 297.0, 313.0, 505.0, 546.0, 441.0};
        assertArrayEquals(conv, DoubleArrays.convolve(a, c), 1e-12);

        conv = new double[]{5.0, 11.0, 16.0, 27.0, 43.0, 70.0, 113.0, 183.0, 126.0};
        assertArrayEquals(conv, DoubleArrays.convolve(a, b), 1e-12);
    }

    @Test
    public void testMax() {
        double[] a = new double[]{-1.0, -0.5, 0.0, 0.5, 1.0};
        assertEquals(1.0, DoubleArrays.max(a), 1e-12);
    }

    @Test
    public void testNorm1() {
        double[] c = {-2.0, 1.0, -0.5, 0.0, 1.5};
        assertEquals(5.0, DoubleArrays.norm1(c), 1e-12);
    }

    @Test
    public void testNorm() {
        double[] c = {-2.0, 1.0, -0.5, 0.0, 1.5};
        assertEquals(Math.sqrt(7.5), DoubleArrays.norm2(c), 1e-12);
    }

    @Test
    public void testNormInf() {
        double[] c = {-2.0, 1.0, -0.5, 0.0, 1.5};
        assertEquals(2.0, DoubleArrays.normInf(c), 1e-12);
    }

    @Test
    public void testNormNegInf() {
        double[] c = {-2.0, 1.0, -0.5, 0.0, 1.5};
        assertEquals(0.0, DoubleArrays.normNegInf(c), 1e-12);
    }

    @Test
    public void testDistance() {
        double[] a = {-1.0, -0.5, 0.0, 0.5, 1.0};
        double[] c = {-2.0, 1.0, -0.5, 0.0, 1.5};
        assertEquals(2, DoubleArrays.distance(a, c), 1e-12);
    }

    @Test
    public void testConcat() {
        double[] a = {-1.0, -0.5, 0.0, 0.5, 1.0};
        double[] c = {-2.0, 1.0, -0.5, 0.0, 1.5};
        double[] concat = new double[]{-1.0, -0.5, 0.0, 0.5, 1.0, -2.0, 1.0, -0.5, 0.0, 1.5};
        assertArrayEquals(concat, DoubleArrays.concatenate(a, c), 1e-12);
    }

    @Test
    public void testAscending() {
        double[] a = new double[]{-1.0, -0.5, 0.0, 0.5, 1.0};
        assertTrue(DoubleArrays.isAscending(a));
        assertFalse(DoubleArrays.isAscending(DoubleArrays.reverse(a)));
    }

    @Test
    public void testRepeat() {
        double[] x = new double[]{1.0, 2.0, 3.0};
        double[] xrep = new double[]{1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
        assertArrayEquals(xrep, DoubleArrays.repeat(x, 3), 1e-12);
    }

    @Test
    public void testCumulativeSum() {
        double[] a = {1, 2, 3, 4, 5, 6};
        double[] sum = DoubleArrays.cumulativeSum(a);
        assertArrayEquals(new double[]{1.0, 3.0, 6.0, 10.0, 15.0, 21.0}, sum, 1e-12);
    }

    @Test
    public void testSum() {
        double[] a = {1, 2, 3, 4, 5, 6};
        double sum = DoubleArrays.sum(a);
        assertEquals(21.0, sum, 1e-12);
    }

    @Test
    public void testSumSquares() {
        double[] a = {1, 2, 3, 4, 5, 6};
        double sum = DoubleArrays.sumSquares(a);
        assertEquals(91.0, sum, 1e-12);
    }

    @Test
    public void testRMS() {
        double[] a = {1, 2, 3, 4, 5, 6};
        double rms = DoubleArrays.rms(a);
        assertEquals(3.894440481849308, rms, 1e-12);
    }

    @Test
    public void testReverse() {
        double[] a = {1, 2, 3, 4, 5};
        double[] expected = {5, 4, 3, 2, 1};

        assertArrayEquals(expected, DoubleArrays.reverse(a), 1e-12);

        a = new double[]{1, 2, 3, 4, 5, 6};
        expected = new double[]{6, 5, 4, 3, 2, 1};

        assertArrayEquals(expected, DoubleArrays.reverse(a), 1e-12);
    }

    @Test
    public void testSubtraction() {
        double[] a = {1, 2, 3, 4, 5};
        double[] b = {10, 9, 8, 7, 6};
        double[] expected = {-9, -7, -5, -3, -1};

        assertArrayEquals(expected, DoubleArrays.subtractElementWise(a, b), 1e-12);

        expected = new double[]{-3, -2, -1, 0, 1};
        assertArrayEquals(expected, DoubleArrays.subtractElementWise(b[3], b), 1e-12);

        expected = new double[]{3, 2, 1, 0, -1};
        assertArrayEquals(expected, DoubleArrays.subtractElementWise(b, b[3]), 1e-12);
    }

    @Test
    public void testAddition() {
        double[] a = {1, 2, 3, 4, 5};
        double[] b = {10, 9, 8, 7, 6};
        double[] expected = {11, 11, 11, 11, 11};

        assertArrayEquals(expected, DoubleArrays.addElementWise(a, b), 1e-12);


        expected = new double[]{17, 16, 15, 14, 13};
        assertArrayEquals(expected, DoubleArrays.addElementWise(b, b[3]), 1e-12);
    }

    @Test
    public void testArgSort() {
        double[] a = {4, 3, 6, -2, 5};
        int[] expected = {3, 1, 0, 4, 2};
        assertArrayEquals(expected, DoubleArrays.argSort(a));
    }

    @Test
    public void testGradient() {
        double[] a = {4, 3, 6, -2, 5};
        double[] expected = {-1, 1, -2.5, -0.5, 7};

        assertArrayEquals(expected, DoubleArrays.gradient(a), 1e-12);

        double[] h = {1, 2, 3, 4, 5};
        assertArrayEquals(expected, DoubleArrays.gradient(a, h), 1e-12);
    }

    @Test
    public void testDotProduct() {
        double[] a = {1, 2, 3, 4, 5};
        double[] b = {10, 9, 8, 7, 6};
        double expected = 110;
        assertEquals(expected, DoubleArrays.dot(a, b), 1e-12);
    }

    @Test
    public void testDivideElementWise() {
        double[] a = {2, 4, 6, 8};
        double d = 2;

        double[] expected = {1, 2, 3, 4};
        double[] actual = DoubleArrays.divideElementWise(a, d);
        assertArrayEquals(expected, actual, 0);
    }
}
