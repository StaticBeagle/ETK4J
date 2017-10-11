package com.wildbitsfoundry.etk4j.math.interpolation;

import static org.junit.Assert.assertEquals;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import com.wildbitsfoundry.etk4j.math.interpolation.CardinalSpline;

public class CardinalSplineTest {
	double[] x;
	double[] y;
	double[] xi;
	double left;
	double right;
	
	@Rule
	public final ExpectedException exception = ExpectedException.none();

	@Before
	public void setArrays() {
		x = new double[] { 0.9, 1.3, 1.9, 2.1 };
		y = new double[] { 1.3, 1.5, 1.85, 2.1 };
		xi = new double[] {1.15, 1.8, 2.0};
		left = -0.5;
		right = 3.0;
	}
	
	@Test
	public void testUniformCatmullRomSplineInterpolation() {
		CardinalSpline cspline = CardinalSpline.newUniformCatmullRomSpline(x, y);
		
		double yi = cspline.evaluateAt(xi[0]);
		assertEquals(1.4140136718749996, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[1]);
		assertEquals(1.795717592592593, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[2]);
		assertEquals(1.9812500000000002, yi, 1e-12);
	}
	
	@Test
	public void testCentripetalCatmullRomSplineInterpolation() {
		CardinalSpline cspline = CardinalSpline.newCentripetalCatmullRomSpline(x, y);
		
		double yi = cspline.evaluateAt(xi[0]);
		assertEquals(1.4235463009194955, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[1]);
		assertEquals(1.7758688274958532, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[2]);
		assertEquals(1.9722494634552241, yi, 1e-12);
	}
	
	@Test
	public void testChordalCatmullRomSplineInterpolation() {
		CardinalSpline cspline = CardinalSpline.newChordalCatmullRomSpline(x, y);
		
		double yi = cspline.evaluateAt(xi[0]);
		assertEquals(1.4273575841257184, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[1]);
		assertEquals(1.71490474060342, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[2]);
		assertEquals(1.9708349756006007, yi, 1e-12);
	}
	
	@Test
	public void testNaturalSplineExtrapolateLeft() {
		CardinalSpline cspline = CardinalSpline.newCentripetalCatmullRomSpline(x, y);

		exception.expect(IndexOutOfBoundsException.class);
		@SuppressWarnings("unused")
		double yi = cspline.evaluateAt(left);
	}
	
	@Test
	public void testNaturalSplineExtrapolateRight() {
		CardinalSpline cspline = CardinalSpline.newCentripetalCatmullRomSpline(x, y);

		exception.expect(IndexOutOfBoundsException.class);
		@SuppressWarnings("unused")
		double yi = cspline.evaluateAt(right);
	}
}
