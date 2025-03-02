package com.wildbitsfoundry.etk4j.statistic.regression;

import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class LinearRegressionTest {
	LinearRegression lr;

	@Before
	public void before() {
		double[] x = { 1, 2, 3, 4, 5.06, 6 };
		double[] y = { 1, 3, 2, 6, 7.23, 7 };
		lr = new LinearRegression(x, y);
	}

	@Test
	public void testBeta() {
		double[] beta = { 1.329893117683651, -0.2962581764029481 };
		assertArrayEquals(beta, lr.beta(), 1e-12);
	}

	@Test
	public void testResiduals() {
		double[] residuals = { -0.03363494128070288, 0.6364719410356461, -1.693421176648005,
				0.976685705668344, 0.7969990009236749, -0.6831005296989581 };
		assertArrayEquals(residuals, lr.residuals(), 1e-12);
	}
	
	@Test
	public void testNormOfResiduals() {
		assertEquals(2.308603870594901, lr.normOfResiduals(), 1e-12);
	}
	
	@Test
	public void testR() {
		assertEquals(0.9243361628117032, lr.R(), 1e-12);
	}

	@Test
	public void testR2() {
		assertEquals(0.8543973418814634, lr.R2(), 1e-12);
	}
	
	@Test
	public void testSSE() {
		assertEquals(5.329651831325758, lr.SSE(), 1e-12);
	}
	
	@Test
	public void testSSR() {
		assertEquals(36.60408333333334, lr.SSR(), 1e-12);
	}
}
