package com.wildbitsfoundry.etk4j.math.interpolation;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import com.wildbitsfoundry.etk4j.math.functions.common.ExtrapolationMethod;
import com.wildbitsfoundry.etk4j.math.interpolation.CubicSpline;

public class CubicSplineTest {
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
	public void testNaturalSplineInterpolation() {
		CubicSpline cspline = CubicSpline.newNaturalSpline(x, y);
		
		double yi = cspline.evaluateAt(xi[0]);
		assertEquals(1.4307218309859153, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[1]);
		assertEquals(1.755721830985916, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[2]);
		assertEquals(1.9682218309859156, yi, 1e-12);
	}
	
	@Test
	public void testParabolicallyTerminatedSplineInterpolation() {
		CubicSpline cspline = CubicSpline.newParabolicallyTerminatedSpline(x, y);
		
		double yi = cspline.evaluateAt(xi[0]);
		assertEquals(1.4321022727272725, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[1]);
		assertEquals(1.7594696969696972, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[2]);
		assertEquals(1.9632575757575759, yi, 1e-12);
	}
	
	@Test
	public void testClampedSplineInterpolation() {
		CubicSpline cspline = CubicSpline.newClampedSpline(x, y, 2, 1);
		
		double yi = cspline.evaluateAt(xi[0]);
		assertEquals(1.5100360576923078, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[1]);
		assertEquals(1.7361111111111118, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[2]);
		assertEquals(1.9814102564102565, yi, 1e-12);
	}
	
	@Test
	public void testNotAKnotSplineInterpolation() {
		CubicSpline cspline = CubicSpline.newNotAKnotSpline(x, y);
		
		double yi = cspline.evaluateAt(xi[0]);
		assertEquals(1.4394531249999998, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[1]);
		assertEquals(1.7593750000000004, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[2]);
		assertEquals(1.9622916666666668, yi, 1e-12);
	}
	
	@Test
	public void testAkimaSplineInterpolation() {
		double[] x = { 0.5, 0.9, 1.3, 1.9, 2.1, 2.2 };
		double[] y = { 1.0, 1.3, 1.5, 1.85, 2.1, 2.4 };
		CubicSpline cspline = CubicSpline.newAkimaSpline(x, y);
		
		double yi = cspline.evaluateAt(xi[0]);
		assertEquals(1.4258655894886363, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[1]);
		assertEquals(1.7887205387205394, yi, 1e-12);
		
		yi = cspline.evaluateAt(xi[2]);
		assertEquals(1.9470219435736678, yi, 1e-12);
	}

	@Test
	public void testNaturalSplineExtrapolateLeft() {
		CubicSpline cspline = CubicSpline.newNaturalSpline(x, y);

		cspline.setExtrapolationMethod(ExtrapolationMethod.ClampToEndPoint);
		double yi = cspline.evaluateAt(left);
		assertEquals("Natural Spline ClampToEndPoint lower bound extrapolation", 1.3, yi, 0.0);

		cspline.setExtrapolationMethod(ExtrapolationMethod.ClampToNaN);
		yi = cspline.evaluateAt(left);
		assertTrue("Natural Spline ClampToEndNaN lower bound extrapolation", Double.isNaN(yi));

		cspline.setExtrapolationMethod(ExtrapolationMethod.ClampToZero);
		yi = cspline.evaluateAt(left);
		assertEquals("Natural Spline ClampToZero lower bound extrapolation", 0.0, yi, 0.0);

		cspline.setExtrapolationMethod(ExtrapolationMethod.Linear);
		yi = cspline.evaluateAt(left);
		assertEquals("Natural Spline Linear lower bound extrapolation", 0.5474178403755874, yi, 1e-12);
		
		cspline.setExtrapolationMethod(ExtrapolationMethod.Natural);
		yi = cspline.evaluateAt(left);
		assertEquals("Natural Spline Natural lower bound extrapolation", 1.1915492957746414, yi, 1e-12);
		
		cspline.setExtrapolationMethod(ExtrapolationMethod.Periodic);
		yi = cspline.evaluateAt(left);
		assertEquals("Natural Spline Periodic lower bound extrapolation", 1.8500000000000005, yi, 1e-12);
		
		cspline.setExtrapolationMethod(ExtrapolationMethod.Throw);
		exception.expect(IndexOutOfBoundsException.class);
		yi = cspline.evaluateAt(left);
	}
	
	@Test
	public void testNaturalSplineExtrapolateRight() {
		CubicSpline cspline = CubicSpline.newNaturalSpline(x, y);

		cspline.setExtrapolationMethod(ExtrapolationMethod.ClampToEndPoint);
		double yi = cspline.evaluateAt(right);
		assertEquals("Natural Spline ClampToEndPoint upper bound extrapolation", 2.1, yi, 0.0);

		cspline.setExtrapolationMethod(ExtrapolationMethod.ClampToNaN);
		yi = cspline.evaluateAt(right);
		assertTrue("Natural Spline ClampToEndNaN upper bound extrapolation", Double.isNaN(yi));

		cspline.setExtrapolationMethod(ExtrapolationMethod.ClampToZero);
		yi = cspline.evaluateAt(right);
		assertEquals("Natural Spline ClampToZero upper bound extrapolation", 0.0, yi, 0.0);

		cspline.setExtrapolationMethod(ExtrapolationMethod.Linear);
		yi = cspline.evaluateAt(right);
		assertEquals("Natural Spline Linear upper bound extrapolation", 3.306338028169013, yi, 1e-12);
		
		cspline.setExtrapolationMethod(ExtrapolationMethod.Natural);
		yi = cspline.evaluateAt(right);
		assertEquals("Natural Spline Natural upper bound extrapolation", 1.6592429577464949, yi, 1e-12);
		
		cspline.setExtrapolationMethod(ExtrapolationMethod.Periodic);
		yi = cspline.evaluateAt(right);
		assertEquals("Natural Spline Periodic upper bound extrapolation", 1.7557218309859157, yi, 1e-12);
		
		cspline.setExtrapolationMethod(ExtrapolationMethod.Throw);
		exception.expect(IndexOutOfBoundsException.class);
		yi = cspline.evaluateAt(right);
	}
}
