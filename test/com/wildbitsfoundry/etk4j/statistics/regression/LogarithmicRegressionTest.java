package com.wildbitsfoundry.etk4j.statistics.regression;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;

public class LogarithmicRegressionTest {
	LogarithmicRegression lr;

	@Before
	public void before() {
		double[] x = { 0.510999, 105.7, 1.78e3 };
		double[] y = { 0, 0.69, 1.06 };
		lr = new LogarithmicRegression(x, y);
	}

	@Test
	public void testBeta() {
		double[] beta = { 0.1298965219335365, 0.08654036624227313 };
		assertArrayEquals(beta, lr.beta(), 1e-12);
	}

	@Test
	public void testResiduals() {
		double[] residuals = { 6.705538069E-4, -1.9367319333130073E-3,  1.266178126E-3 };
		assertArrayEquals(residuals, lr.residuals(), 1e-12);
	}

	@Test
	public void testNormOfResiduals() {
		assertEquals(0.0024091035754607677, lr.normOfResiduals(), 1e-12);
	}

	@Test
	public void testR2() {
		assertEquals(0.9999899738914477, lr.R2(), 1e-12);
	}

	@Test
	public void testSSE() {
		assertEquals(5.803780037297855E-6, lr.SSE(), 1e-12);
	}

	@Test
	public void testSSR() {
		assertEquals(0.57886666666666668, lr.SSR(), 1e-12);
	}
}
