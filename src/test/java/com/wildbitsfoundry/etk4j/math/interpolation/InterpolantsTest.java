package com.wildbitsfoundry.etk4j.math.interpolation;

import org.junit.Test;

import static com.wildbitsfoundry.etk4j.math.interpolation.Interpolation.*;
import static org.junit.Assert.assertEquals;

public class InterpolantsTest {

	@Test
	public void testLinear() {
		double x = linear(1, 2, 1, 2, 1.5);
		assertEquals(1.5, x, 1e-12);
	}
	
	@Test
	public void testNeville() {
		double x = neville(new double[] {1, 2}, new double[] {1, 2}, 1.5);
		assertEquals(1.5, x, 1e-12);
	}
	
	@Test
	public void testParabola() {
		double x = quadratic(1, 2, 3, 1, 4, 9, 2.5);
		assertEquals(6.25, x, 1e-12);
	}
	
	@Test
	public void testSpline() {
		double[] x = {1, 2};
		double[] y = {1, 2};
		
		double d = spline(x, y, 1.5);
		assertEquals(1.5, d, 1e-12);
		
		x = new double[] {1, 2, 3};
		y = new double[] {1, 4, 9};
		
		d = spline(x, y, 2.5);
		assertEquals(6.25, d, 1e-12);
		
		x = new double[] {1, 2, 3, 4};
		y = new double[] {1, 8, 27, 64};
		
		d = spline(x, y, 3.5);
		assertEquals(42.875, d, 1e-12);
	}
}
