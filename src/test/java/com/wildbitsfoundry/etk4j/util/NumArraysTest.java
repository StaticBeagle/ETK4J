package com.wildbitsfoundry.etk4j.util;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Before;
import org.junit.Test;

public class NumArraysTest {
	double[] a;
	double[] b;
	double[] c;

	@Before
	public void before() {
		a = new double[] { -1.0, -0.5, 0.0, 0.5, 1.0 };
		b = new double[] { 0.1000, 0.31622776601683794, 1.0000, 3.1622776601683795, 10.0000 };
		c = new double[] { -2.0, 1.0, -0.5, 0.0, 1.5 };
	}

	@Test
	public void testLinspace() {
		assertArrayEquals(a, NumArrays.linSpace(-1.0, 1.0, 5), 1e-12);
	}

	@Test
	public void testLinspaceStep() {
		assertArrayEquals(a, NumArrays.linSteps(-1.0, 1.0, 0.5), 1e-12);
	}

	@Test
	public void testLogspace() {
		assertArrayEquals(b, NumArrays.logSpace(-1, 1, 5), 1e-12);
	}

	@Test
	public void testConv() {
		double[] conv = new double[] { 1.0000, 1.0000, 0.2500, -1.0000, -2.5000, -1.0000, 0.2500, 1.0000, 1.0000 };
		assertArrayEquals(conv, NumArrays.convolution(a, a), 1e-12);
	}

	@Test
	public void testMax() {
		assertEquals(1.0, NumArrays.max(a), 1e-12);
	}

	@Test
	public void testNorm1() {
		assertEquals(5.0, NumArrays.norm1(c), 1e-12);
	}

	@Test
	public void testNorm() {
		assertEquals(Math.sqrt(7.5), NumArrays.norm2(c), 1e-12);
	}

	@Test
	public void testNormInf() {
		assertEquals(2.0, NumArrays.normInf(c), 1e-12);
	}

	@Test
	public void testConcat() {
		double[] concat = new double[] { -1.0, -0.5, 0.0, 0.5, 1.0, -2.0, 1.0, -0.5, 0.0, 1.5 };
		assertArrayEquals(concat, NumArrays.concatenate(a, c), 1e-12);
	}

	@Test
	public void testAscending() {
		assertTrue(NumArrays.isAscending(a));
		assertFalse(NumArrays.isAscending(NumArrays.reverse(a)));
	}

	@Test
	public void testRepeat() {
		double[] x = new double[] { 1.0, 2.0, 3.0 };
		double[] xrep = new double[] { 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
		assertArrayEquals(xrep, NumArrays.repeat(x, 3), 1e-12);
	}
	
	@Test
	public void testCummulativeSum() {
		double[] a = {1, 2, 3, 4 ,5 ,6};
		double[] sum = NumArrays.cumulativeSum(a);
		assertArrayEquals(new double[] {1.0, 3.0, 6.0, 10.0, 15.0, 21.0}, sum, 1e-12);
	}
	
	@Test
	public void testSum() {
		double[] a = {1, 2, 3, 4 ,5 ,6};
		double sum = NumArrays.sum(a);
		assertEquals(21.0, sum, 1e-12);
	}
	
	@Test
	public void testRMS() {
		double[] a = {1, 2, 3, 4 ,5 ,6};
		double rms = NumArrays.rms(a);
		assertEquals(3.894440481849308, rms, 1e-12);
	}
}
