package com.wildbitsfoundry.etk4j.math.polynomials;

import static org.junit.Assert.assertArrayEquals;

import org.junit.Before;
import org.junit.Test;

public class PolynomialsTest {
	double[] a;
	double[] b;
	double[] c;

	@Before
	public void before() {
		a = new double[] { -1.0, -0.5, 0.0, 0.5, 1.0 };
		b = new double[] { 0.1000, 0.31622776601683794, 1.0000, 3.1622776601683795, 10.0000 };
		c = new double[] { -2.0, 1.0, -0.5, 0.0, 1.5 };
	}

	@Test
	public void testPolyfit() {
		double[] x = new double[] { 1, 2, 3, 4 };
		double[] y = new double[] { 1, 10, 12, 15 };

		// Over determined
		Polynomial poly = Polynomial.polyfit(x, y, 2);
		assertArrayEquals(new double[] { -1.5, 11.9, -9.0 }, poly.getCoefficients(), 1e-12);

		// Unique solution
		Polynomial poly2 = Polynomial.polyfit(x, y, 3);
		assertArrayEquals(new double[] { 1.333333333333, -11.5, 34.1666666666664, -23.0 }, poly2.getCoefficients(),
				1e-12);

		// Under determined
		Polynomial poly3 = Polynomial.polyfit(x, y, 4);
		assertArrayEquals(new double[] { 0.6079433590792012, -4.746100257458709, 9.778017567772217, 3.769498712706295,
				-8.40935938209905 }, poly3.getCoefficients(), 1e-12);
	}
}
