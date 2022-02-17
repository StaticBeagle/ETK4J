package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.extrapolation.ExtrapolationMethod;
import com.wildbitsfoundry.etk4j.util.NumArrays;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import static org.junit.Assert.*;

public class LinearSplineTest {
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
	public void testLinearSplineInterpolation() {
		LinearSpline ls = LinearSpline.newLinearSpline(x, y);
		
		double yi = ls.evaluateAt(xi[0]);
		assertEquals(1.425, yi, 1e-12);
		
		yi = ls.evaluateAt(xi[1]);
		assertEquals(1.7916666666666667, yi, 1e-12);
		
		yi = ls.evaluateAt(xi[2]);
		assertEquals(1.975, yi, 1e-12);

		double[] expected = {0.0, 22.704, 45.408, 68.112, 90.816, 113.52000000000001, 136.224, 158.928, 181.632,
				204.336, 227.04, 254.188, 281.336, 308.484, 335.63199999999995, 362.78, 393.69399999999996, 424.608,
				455.522, 486.43600000000004, 517.35, 551.5980000000001, 585.846, 622.8833333333333, 662.71,
				702.5366666666666, 742.3633333333333, 782.19, 822.0166666666667, 861.8433333333332, 901.67};

		double[] x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
		double[] y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
		ls = LinearSpline.newLinearSpline(x, y);
		assertArrayEquals(expected, ls.evaluateAt(NumArrays.linSteps(0, 30)), 1e-12);

		assertEquals(393.69399999999996, ls.evaluateAt(16.0), 1e-12);
		assertEquals(30.91400000000001, ls.differentiate(16.0), 1e-12);
		assertEquals(30.91400000000001, ls.evaluateDerivativeAt(2, 16), 1e-12);
		assertEquals(102.16799999999999, ls.evaluateAntiDerivativeAt(0, 3.0), 1e-12);
		assertEquals(102.16799999999999, ls.integrate(3.0), 1e-12);
		assertEquals(102.16799999999999, ls.integrate(0.0, 3.0), 1e-12);
		assertEquals(1612.173, ls.integrate(11.,16.0), 1e-12);
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
		
		lspline.setExtrapolationMethod(ExtrapolationMethod.THROW);
		exception.expect(IndexOutOfBoundsException.class);
		yi = lspline.evaluateAt(right);
	}
}
