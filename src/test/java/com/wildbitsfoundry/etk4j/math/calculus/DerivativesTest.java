package com.wildbitsfoundry.etk4j.math.calculus;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import com.wildbitsfoundry.etk4j.math.function.UnivariateFunction;

public final class DerivativesTest {
	
	@Test
	public void testForwardDifference() {
		UnivariateFunction func = x -> x * x + 2 * x + 1;
		double x = 1.4;
		double h = 1e-7;
		double fp = Derivatives.forwardDifference(func, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testForwardDifference3Points() {
		UnivariateFunction func = x -> x * x + 2 * x + 1;
		double x = 1.4;
		double h = 1e-7;
		double fp = Derivatives.forwardDifference3Points(func, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testForwardDifference5Points() {
		UnivariateFunction func = x -> x * x + 2 * x + 1;
		double x = 1.4;
		double h = 1e-7;
		double fp = Derivatives.forwardDifference5Points(func, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testBackwardDifference() {
		UnivariateFunction func = x -> x * x + 2 * x + 1;
		double x = 1.4;
		double h = 1e-7;
		double fp = Derivatives.backwardDifference(func, x, h);
		assertEquals(4.8, fp, 1e-6);
	}
	
	@Test
	public void testBackwardDifference3Points() {
		UnivariateFunction func = x -> x * x + 2 * x + 1;
		double x = 1.4;
		double h = 1e-7;
		double fp = Derivatives.backwardDifference3Points(func, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testBackwardDifference5Points() {
		UnivariateFunction func = x -> x * x + 2 * x + 1;
		double x = 1.4;
		double h = 1e-7;
		double fp = Derivatives.backwardDifference5Points(func, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testCenteredDifference() {
		UnivariateFunction func = x -> x * x + 2 * x + 1;
		double x = 1.4;
		double h = 1e-7;
		double fp = Derivatives.centeredDifference(func, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testCenteredDifference5Points() {
		UnivariateFunction func = x -> x * x + 2 * x + 1;
		double x = 1.4;
		double h = 1e-7;
		double fp = Derivatives.centeredDifference5Points(func, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testCenteredDifference7Points() {
		UnivariateFunction func = x -> x * x + 2 * x + 1;
		double x = 1.4;
		double h = 1e-7;
		double fp = Derivatives.centeredDifference7Points(func, x, h);
		assertEquals(4.8, fp, 1e-7);
	}

	@Test
	public void testRidders() {
		UnivariateFunction func = x -> x * x + 2 * x + 1;
		double x = 1.4;
		double h = 1e-7;
		double fp = Derivatives.ridders(func, x, h);
		assertEquals(4.8, fp, 1e-7);
	}

}