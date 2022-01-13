package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.util.NumArrays;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import java.util.Arrays;

import static org.junit.Assert.*;

public class CubicSplineTest {
	static double[] x;
	static double[] y;
	static double[] xi;
	static double left;
	static double right;
//
//	// TODO test derivatives and definite/indefinite integrals
//	// add tests for quadratic spline
//
	@Rule
	public final ExpectedException exception = ExpectedException.none();

	@BeforeClass
	public static void setArrays() {
		x = new double[] { 0.9, 1.3, 1.9, 2.1 };
		y = new double[] { 1.3, 1.5, 1.85, 2.1 };
		xi = new double[] {1.15, 1.8, 2.0};
		left = -0.5;
		right = 3.0;
	}

	@Test
	public void testNaturalSplineInterpolation() {
		CubicSpline cs = CubicSpline.newNaturalSpline(x, y);

		double yi = cs.evaluateAt(xi[0]);
		assertEquals(1.4307218309859153, yi, 1e-12);

		yi = cs.evaluateAt(xi[1]);
		assertEquals(1.755721830985916, yi, 1e-12);

		yi = cs.evaluateAt(xi[2]);
		assertEquals(1.9682218309859156, yi, 1e-12);

		double[] expected = {0.0, 21.438503269035536, 42.95370330964467, 64.62229689340103, 86.52098079187817,
				108.72645177664975, 131.31540661928935, 154.36454209137057, 177.95055496446702, 202.15014201015228,
				227.04, 252.68284377664975, 279.0854604670051, 306.2406552690355, 334.1412333807106, 362.78,
				392.15420158375633, 422.2788496243655, 453.1733968730964, 484.8572960812183, 517.35, 550.7171732791878,
				585.2093281624366, 621.1126457394247, 658.4708168934011, 697.0850423011844, 736.7459795871405,
				777.2442863756345, 818.3706202910322, 859.9156389576988, 901.67};

		double[] x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
		double[] y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
		cs = CubicSpline.newNaturalSpline(x, y);
		assertArrayEquals(expected, cs.evaluateAt(NumArrays.linSteps(0, 30)), 1e-12);

		assertEquals(392.15420158375633, cs.evaluateAt(16.0), 1e-12);
		assertEquals(29.746182686971242, cs.differentiate(16.0), 1e-12);
		assertEquals(29.746182686971242, cs.evaluateDerivativeAt(2, 16), 1e-12);
		assertEquals(377.4053741184433, cs.evaluateAntiDerivativeAt(2, 16.0), 1e-12);
		assertEquals(2947.3965035600677, cs.integrate(16.0), 1e-12);
		assertEquals(2947.3965035600677, cs.integrate(0.0, 16.0), 1e-12);
		assertEquals(1604.3556840203046, cs.integrate(11.,16.0), 1e-12);
	}
//
//	@Test
//	public void testParabolicallyTerminatedSplineInterpolation() {
//		CubicSpline cspline = CubicSpline.newParabolicallyTerminatedSpline(x, y);
//
//		double yi = cspline.evaluateAt(xi[0]);
//		assertEquals(1.4321022727272725, yi, 1e-12);
//
//		yi = cspline.evaluateAt(xi[1]);
//		assertEquals(1.7594696969696972, yi, 1e-12);
//
//		yi = cspline.evaluateAt(xi[2]);
//		assertEquals(1.9632575757575759, yi, 1e-12);
//	}
//
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
		CubicSpline cs = CubicSpline.newNotAKnotSpline(x, y);

		double yi = cs.evaluateAt(xi[0]);
		assertEquals(1.4394531249999998, yi, 1e-12);

		yi = cs.evaluateAt(xi[1]);
		assertEquals(1.7593750000000004, yi, 1e-12);

		yi = cs.evaluateAt(xi[2]);
		assertEquals(1.9622916666666668, yi, 1e-12);

		double[] expected = {0.0, 20.514440000000018, 41.45489777777781, 62.844080000000034, 84.70469333333337,
				107.05944444444448, 129.93104000000002, 153.3421866666667, 177.31559111111113, 201.87396, 227.04,
				252.83641777777777, 279.28592, 306.4112133333333, 334.2350044444444, 362.78, 392.0707644444444,
				422.13929333333334, 453.01944000000003, 484.7450577777778, 517.35, 550.87058, 585.3529511111111,
				620.8457266666667, 657.39752, 695.0569444444445, 733.8726133333333, 773.89314, 815.1671377777777,
				857.7432200000001, 901.67};

		double[] x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
		double[] y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
		cs = CubicSpline.newNotAKnotSpline(x, y);
		assertArrayEquals(expected, cs.evaluateAt(NumArrays.linSteps(0, 30)), 1e-12);

		assertEquals(392.0707644444444, cs.evaluateAt(16.0), 1e-12);
		assertEquals(29.674004444444453, cs.differentiate(16.0), 1e-12);
		assertEquals(29.674004444444453, cs.evaluateDerivativeAt(2, 16), 1e-12);
		assertEquals(377.36197907407404, cs.evaluateAntiDerivativeAt(2, 16.0), 1e-12);
		assertEquals(2936.883854074074, cs.integrate(16.0), 1e-12);
		assertEquals(2936.883854074074, cs.integrate(0.0, 16.0), 1e-12);
		assertEquals(1604.869493148148, cs.integrate(11.,16.0), 1e-12);
	}

//	@Test
//	public void testAkimaSplineInterpolation() {
//		double[] x = { 0.5, 0.9, 1.3, 1.9, 2.1, 2.2 };
//		double[] y = { 1.0, 1.3, 1.5, 1.85, 2.1, 2.4 };
//		CubicSpline cspline = CubicSpline.newAkimaSpline(x, y);
//
//		double yi = cspline.evaluateAt(xi[0]);
//		assertEquals(1.4258655894886363, yi, 1e-12);
//
//		yi = cspline.evaluateAt(xi[1]);
//		assertEquals(1.7887205387205394, yi, 1e-12);
//
//		yi = cspline.evaluateAt(xi[2]);
//		assertEquals(1.9470219435736678, yi, 1e-12);
//	}
//
	@Test
	public void testNaturalSplineExtrapolateLeft() {
		CubicSpline cspline = CubicSpline.newNaturalSpline(x, y);

		cspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_END_POINT);
		double yi = cspline.evaluateAt(left);
		assertEquals("Natural Spline ClampToEndPoint lower bound extrapolation", 1.3, yi, 0.0);

		cspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_NAN);
		yi = cspline.evaluateAt(left);
		assertTrue("Natural Spline ClampToEndNaN lower bound extrapolation", Double.isNaN(yi));

		cspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_ZERO);
		yi = cspline.evaluateAt(left);
		assertEquals("Natural Spline ClampToZero lower bound extrapolation", 0.0, yi, 0.0);

		cspline.setExtrapolationMethod(ExtrapolationMethod.LINEAR);
		yi = cspline.evaluateAt(left);
		assertEquals("Natural Spline Linear lower bound extrapolation", 0.5474178403755874, yi, 1e-12);

		cspline.setExtrapolationMethod(ExtrapolationMethod.NATURAL);
		yi = cspline.evaluateAt(left);
		assertEquals("Natural Spline Natural lower bound extrapolation", 1.1915492957746414, yi, 1e-12);

		cspline.setExtrapolationMethod(ExtrapolationMethod.PERIODIC);
		yi = cspline.evaluateAt(left);
		assertEquals("Natural Spline Periodic lower bound extrapolation", 1.8500000000000005, yi, 1e-12);

		cspline.setExtrapolationMethod(ExtrapolationMethod.THROW);
		exception.expect(IndexOutOfBoundsException.class);
		cspline.evaluateAt(left);
	}

	@Test
	public void testNaturalSplineExtrapolateRight() {
		CubicSpline cspline = CubicSpline.newNaturalSpline(x, y);

		cspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_END_POINT);
		double yi = cspline.evaluateAt(right);
		assertEquals("Natural Spline ClampToEndPoint upper bound extrapolation", 2.1, yi, 0.0);

		cspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_NAN);
		yi = cspline.evaluateAt(right);
		assertTrue("Natural Spline ClampToEndNaN upper bound extrapolation", Double.isNaN(yi));

		cspline.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_ZERO);
		yi = cspline.evaluateAt(right);
		assertEquals("Natural Spline ClampToZero upper bound extrapolation", 0.0, yi, 0.0);

		cspline.setExtrapolationMethod(ExtrapolationMethod.LINEAR);
		yi = cspline.evaluateAt(right);
		assertEquals("Natural Spline Linear upper bound extrapolation", 3.306338028169013, yi, 1e-12);

		cspline.setExtrapolationMethod(ExtrapolationMethod.NATURAL);
		yi = cspline.evaluateAt(right);
		assertEquals("Natural Spline Natural upper bound extrapolation", 1.6592429577464949, yi, 1e-12);

		cspline.setExtrapolationMethod(ExtrapolationMethod.PERIODIC);
		yi = cspline.evaluateAt(right);
		assertEquals("Natural Spline Periodic upper bound extrapolation", 1.7557218309859157, yi, 1e-12);

		cspline.setExtrapolationMethod(ExtrapolationMethod.THROW);
		exception.expect(IndexOutOfBoundsException.class);
		cspline.evaluateAt(right);
	}
}
