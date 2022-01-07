package com.wildbitsfoundry.etk4j.math.interpolation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class LinearSplineTest {
	double[] x;
	double[] y;
	double[] xi;
	double left;
	double right;

	// TODO test derivatives and definite/indefinite integrals
	
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
	public void testLinearSplineInterpolation() {
		LinearSpline lspline = LinearSpline.newLinearSpline(x, y);
		
		double yi = lspline.evaluateAt(xi[0]);
		assertEquals(1.425, yi, 1e-12);
		
		yi = lspline.evaluateAt(xi[1]);
		assertEquals(1.7916666666666667, yi, 1e-12);
		
		yi = lspline.evaluateAt(xi[2]);
		assertEquals(1.975, yi, 1e-12);
	}
	
	

	@Test
	public void testLinearSplineExtrapolateLeft() {
		LinearSpline lspline = LinearSpline.newLinearSpline(x, y);

		lspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_END_POINT);
		double yi = lspline.evaluateAt(left);
		assertEquals("Linear Spline ClampToEndPoint lower bound extrapolation", 1.3, yi, 0.0);

		lspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_NAN);
		yi = lspline.evaluateAt(left);
		assertTrue("Linear Spline ClampToEndNaN lower bound extrapolation", Double.isNaN(yi));

		lspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_ZERO);
		yi = lspline.evaluateAt(left);
		assertEquals("Linear Spline ClampToZero lower bound extrapolation", 0.0, yi, 0.0);

		lspline.setExtrapolationMethod(ExtrapolationMethod.LINEAR);
		yi = lspline.evaluateAt(left);
		assertEquals("Linear Spline Linear lower bound extrapolation", 0.6000000000000002, yi, 1e-12);
		
		lspline.setExtrapolationMethod(ExtrapolationMethod.NATURAL);
		yi = lspline.evaluateAt(left);
		assertEquals("Linear Spline Natural lower bound extrapolation", 0.6000000000000002, yi, 1e-12);
		
		lspline.setExtrapolationMethod(ExtrapolationMethod.PERIODIC);
		yi = lspline.evaluateAt(left);
		assertEquals("Linear Spline Periodic lower bound extrapolation", 1.8500000000000005, yi, 1e-12);
		
		lspline.setExtrapolationMethod(ExtrapolationMethod.THROW);
		exception.expect(IndexOutOfBoundsException.class);
		yi = lspline.evaluateAt(left);
	}
	
	@Test
	public void testLinearSplineExtrapolateRight() {
		LinearSpline lspline = LinearSpline.newLinearSpline(x, y);

		lspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_END_POINT);
		double yi = lspline.evaluateAt(right);
		assertEquals("Linear Spline ClampToEndPoint upper bound extrapolation", 2.1, yi, 0.0);

		lspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_NAN);
		yi = lspline.evaluateAt(right);
		assertTrue("Linear Spline ClampToEndNaN upper bound extrapolation", Double.isNaN(yi));

		lspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_ZERO);
		yi = lspline.evaluateAt(right);
		assertEquals("Linear Spline ClampToZero upper bound extrapolation", 0.0, yi, 0.0);

		lspline.setExtrapolationMethod(ExtrapolationMethod.LINEAR);
		yi = lspline.evaluateAt(right);
		assertEquals("Linear Spline Linear upper bound extrapolation", 3.2249999999999988, yi, 1e-12);
		
		lspline.setExtrapolationMethod(ExtrapolationMethod.NATURAL);
		yi = lspline.evaluateAt(right);
		assertEquals("Linear Spline Natural upper bound extrapolation", 3.2249999999999988, yi, 1e-12);
		
		lspline.setExtrapolationMethod(ExtrapolationMethod.PERIODIC);
		yi = lspline.evaluateAt(right);
		assertEquals("Linear Spline Periodic upper bound extrapolation", 1.7916666666666667, yi, 1e-12);
		
		lspline.setExtrapolationMethod(ExtrapolationMethod.THROW);
		exception.expect(IndexOutOfBoundsException.class);
		yi = lspline.evaluateAt(right);
	}
}
