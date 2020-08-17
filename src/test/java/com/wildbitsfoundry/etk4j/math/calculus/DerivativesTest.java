package com.wildbitsfoundry.etk4j.math.calculus;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

public final class DerivativesTest {
	
	UnivariateFunction fn = x -> x * x + 2 * x + 1;
	double x = 1.4;
	double h = 1e-7;
	
	@Test
	public void testForwardDifference() {
		double fp = Derivatives.forwardDifference(fn, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testForwardDifference3Points() {
		double fp = Derivatives.forwardDifference3Points(fn, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testForwardDifference5Points() {
		double fp = Derivatives.forwardDifference5Points(fn, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testBackwardDifference() {
		double fp = Derivatives.backwardDifference(fn, x, h);
		assertEquals(4.8, fp, 1e-6);
	}
	
	@Test
	public void testBackwardDifference3Points() {
		double fp = Derivatives.backwardDifference3Points(fn, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testBackwardDifference5Points() {
		double fp = Derivatives.backwardDifference5Points(fn, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testCenteredDifference() {
		double fp = Derivatives.centeredDifference(fn, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testCenteredDifference5Points() {
		double fp = Derivatives.centeredDifference5Points(fn, x, h);
		assertEquals(4.8, fp, 1e-7);
	}
	
	@Test
	public void testCenteredDifference7Points() {
		double fp = Derivatives.centeredDifference7Points(fn, x, h);
		assertEquals(4.8, fp, 1e-7);
	}

}