package com.wildbitsfoundry.etk4j.math.polynomials;

import static org.junit.Assert.assertArrayEquals;

import org.junit.Before;
import org.junit.Test;

import com.wildbitsfoundry.etk4j.util.ArrayUtils;

public class PolynomialTest {
	double[] a;
	double[] b;
	double[] c;
	
	@Before
	public void before() {
		a = new double[] { -1.0, -0.5, 0.0, 0.5, 1.0 };
		b = new double[] { 0.1000, 0.31622776601683794 , 1.0000, 3.1622776601683795, 10.0000 };
		c = new double[] { -2.0, 1.0, -0.5, 0.0, 1.5 };
	}
	
	@Test
	public void testPolyfit() {
		assertArrayEquals(a, ArrayUtils.linspace(-1.0, 1.0, 5), 1e-12);
	}
}
