package com.wildbitsfoundry.etk4j.math.curvefitting;

import static org.junit.Assert.assertArrayEquals;

import org.junit.Test;

public final class CurveFittingTest {
	
	@Test
	public void testLine() {
		double[] line = CurveFitting.line(0.0, 1.0, 1.0, 2.0);
		assertArrayEquals(new double[] {1.0, 1.0}, line, 1e-12);
	}
	
	@Test
	public void testLinear() {
		double[] x = { 1, 2, 3, 4, 5.06, 6 };
		double[] y = { 1, 3, 2, 6, 7.23, 7 };
		double[] beta = { 1.329893117683651, -0.2962581764029481 };
		double[] linear = CurveFitting.linear(x, y);
		assertArrayEquals(beta, linear, 1e-12);
	}
	
	@Test
	public void testParabola() {
		double[] x = { -4.62, -1.96, 1.89 };
		double[] y = { 2.99, 2.12, 4.02 };
		double[] coefficients = { 0.12604825847610224, 0.5023298715998206, 2.6203395585738543 };
		double[] parabola = CurveFitting.parabola(x[0], x[1], x[2],y[0], y[1], y[2]);
		assertArrayEquals(coefficients, parabola, 1e-12);
	}

	@Test
	public void testPolynomial() {
		double[] x = { -4.62, -1.96, 1.89 };
		double[] y = { 2.99, 2.12, 4.02 };
		double[] coefficients = { 0.12604825847610227, 0.5023298715998206, 2.6203395585738543 };
		double[] quadratic = CurveFitting.polynomial(x, y, 2);
		assertArrayEquals(coefficients, quadratic, 1e-12);
	}
	
	@Test
	public void testExponential() {
		double[] x = { 1, 5, 8 };
		double[] y = { 2, 6, 10 };
		double[] coefficients = { 1.8380687899821915, 0.21612342866514117 };
		double[] exponential = CurveFitting.exponential(x, y);
		assertArrayEquals(coefficients, exponential, 1e-12);
	}
	
	@Test
	public void testLogarithmic() {
		double[] x = { 1, 5, 8 };
		double[] y = { 2, 6, 10 };
		double[] coefficients = { 1.6997797260492338, 3.4971760347075427 };
		double[] log = CurveFitting.logarithmic(x, y);
		assertArrayEquals(coefficients, log, 1e-12);
	}
	
	@Test
	public void testPowerLaw() {
		double[] x = { 1, 5, 8 };
		double[] y = { 2, 6, 10 };
		double[] coefficients = { 1.960117903982167, 0.7504929064884995 };
		double[] pow = CurveFitting.power(x, y);
		assertArrayEquals(coefficients, pow, 1e-12);
	}
	
}
