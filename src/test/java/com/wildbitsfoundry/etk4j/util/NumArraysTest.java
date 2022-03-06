package com.wildbitsfoundry.etk4j.util;

import static com.wildbitsfoundry.etk4j.util.NumArrays.convolution;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;

public class NumArraysTest {


    @Test
    public void testLinspace() {
        double[] a = new double[]{-1.0, -0.5, 0.0, 0.5, 1.0};
        assertArrayEquals(a, NumArrays.linSpace(-1.0, 1.0, 5), 1e-12);
    }

    @Test
    public void testLinspaceStep() {
        double[] a = new double[]{-1.0, -0.5, 0.0, 0.5, 1.0};
        assertArrayEquals(a, NumArrays.linSteps(-1.0, 1.0, 0.5), 1e-12);
    }

    @Test
    public void testLogspace() {
        double[] b = new double[]{0.1000, 0.31622776601683794, 1.0000, 3.1622776601683795, 10.0000};
        assertArrayEquals(b, NumArrays.logSpace(-1, 1, 5), 1e-12);
    }

    @Test
    public void testConv() {
        double[] a = new double[]{1, 1, 2, 3, 5, 8, 13, 21};
        double[] b = new double[]{5, 6};
        double[] c = new double[]{0, 1, 1, 0, 2, 3, 5, 0, 8, 13, 21};

        double[] conv = {1.0, 2.0, 5.0, 10.0, 20.0, 38.0, 71.0, 130.0, 167.0, 242.0, 320.0, 418.0, 505.0, 546.0, 441.0};
        assertArrayEquals(conv, convolution(a, a), 1e-12);

        conv = new double[]{0.0, 1.0, 2.0, 3.0, 7.0, 13.0, 25.0, 38.0, 71.0, 88.0, 125.0, 192.0,
                249.0, 297.0, 313.0, 505.0, 546.0, 441.0};
        assertArrayEquals(conv, convolution(a, c), 1e-12);

        conv = new double[]{5.0, 11.0, 16.0, 27.0, 43.0, 70.0, 113.0, 183.0, 126.0};
        assertArrayEquals(conv, convolution(a, b), 1e-12);
    }

    @Test
    public void testMax() {
        double[] a = new double[]{-1.0, -0.5, 0.0, 0.5, 1.0};
        assertEquals(1.0, NumArrays.max(a), 1e-12);
    }

    @Test
    public void testNorm1() {
        double[] c = {-2.0, 1.0, -0.5, 0.0, 1.5};
        assertEquals(5.0, NumArrays.norm1(c), 1e-12);
    }

    @Test
    public void testNorm() {
        double[] c = {-2.0, 1.0, -0.5, 0.0, 1.5};
        assertEquals(Math.sqrt(7.5), NumArrays.norm2(c), 1e-12);
    }

    @Test
    public void testNormInf() {
        double[] c = {-2.0, 1.0, -0.5, 0.0, 1.5};
        assertEquals(2.0, NumArrays.normInf(c), 1e-12);
    }

    @Test
    public void testConcat() {
        double[] a = {-1.0, -0.5, 0.0, 0.5, 1.0};
        double[] c = {-2.0, 1.0, -0.5, 0.0, 1.5};
        double[] concat = new double[]{-1.0, -0.5, 0.0, 0.5, 1.0, -2.0, 1.0, -0.5, 0.0, 1.5};
        assertArrayEquals(concat, NumArrays.concatenate(a, c), 1e-12);
    }

    @Test
    public void testAscending() {
        double[] a = new double[]{-1.0, -0.5, 0.0, 0.5, 1.0};
        assertTrue(NumArrays.isAscending(a));
        assertFalse(NumArrays.isAscending(NumArrays.reverse(a)));
    }

    @Test
    public void testRepeat() {
        double[] x = new double[]{1.0, 2.0, 3.0};
        double[] xrep = new double[]{1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
        assertArrayEquals(xrep, NumArrays.repeat(x, 3), 1e-12);
    }

    @Test
    public void testCummulativeSum() {
        double[] a = {1, 2, 3, 4, 5, 6};
        double[] sum = NumArrays.cumulativeSum(a);
        assertArrayEquals(new double[]{1.0, 3.0, 6.0, 10.0, 15.0, 21.0}, sum, 1e-12);
    }

    @Test
    public void testSum() {
        double[] a = {1, 2, 3, 4, 5, 6};
        double sum = NumArrays.sum(a);
        assertEquals(21.0, sum, 1e-12);
    }

    @Test
    public void testRMS() {
        double[] a = {1, 2, 3, 4, 5, 6};
        double rms = NumArrays.rms(a);
        assertEquals(3.894440481849308, rms, 1e-12);
    }

    @Test
    public void testReverse() {
        double[] a = {1, 2, 3, 4, 5};
        double[] expected = {5, 4, 3, 2, 1};

        assertArrayEquals(expected, NumArrays.reverse(a), 1e-12);

        a = new double[]{1, 2, 3, 4, 5, 6};
        expected = new double[]{6, 5, 4, 3, 2, 1};

        assertArrayEquals(expected, NumArrays.reverse(a), 1e-12);
    }
}
