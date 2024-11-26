package com.wildbitsfoundry.etk4j.math.calculus;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

public final class IntegralsTest {
	
	
	@Test
	public void testCummulativeTrapz() {
		double[] a = {1, 2, 3, 4 ,5 ,6};
		double[] trapz = Integrals.cumulativeTrapz(a);
		assertArrayEquals(new double[] {0.0, 1.5, 4.0, 7.5, 12.0, 17.5}, trapz, 1e-12);
	}
	
	@Test
	public void testTrapz() {
		double[] a = {1, 2, 3, 4 ,5 ,6};
		double trapz = Integrals.trapz(a);
		assertEquals(17.5, trapz, 1e-12);
		
		UnivariateFunction fx = x -> Math.sin(x * x);
		trapz = Integrals.trapz(fx, 0, Math.PI / 2.0, 1000);
		assertEquals(0.8281158242079556, trapz, 1e-12);
	}	
	
	@Test
	public void testSimpson() {
		UnivariateFunction fx = x -> Math.sin(x * x);
		double simp = Integrals.simpson(fx, 0, Math.PI / 2.0, 1000);
		assertEquals(0.8281163288433171, simp, 1e-12);
	}

	@Test
	public void testGaussianQuadrature() {
		double[] e = new double[2];
		UnivariateFunction fx = x -> Math.sin(x * x);
		e[0] = ConstantsETK.FLOAT_EPS; // relative tol
		e[1] = ConstantsETK.FLOAT_EPS; // absolute tol
		double qadrat = Integrals.gaussianQuadrature(fx, 0, Math.PI / 2.0);
		assertEquals(0.8281262542785998, qadrat, 1e-12);

		fx = x -> Math.sin(x * x) / x;
		qadrat = Integrals.gaussianQuadrature(fx, Math.PI / 4.0, Math.PI / 2.0);
		assertEquals(0.5832680904501192, qadrat, 1e-12);
	}
	
	@Test
	public void testAdaptiveGaussianQuadrature() {
		double[] e = new double[2];
		UnivariateFunction fx = x -> Math.sin(x * x);
		e[0] = ConstantsETK.FLOAT_EPS; // relative tol
		e[1] = ConstantsETK.FLOAT_EPS; // absolute tol
		double qadrat = Integrals.adaptiveGaussianQuadrature(fx, 0, Math.PI / 2.0, e[0], e[1], 150);
		assertEquals(0.8281163273860381, qadrat, 1e-12);

		fx = x -> Math.sin(x * x) / x;
		qadrat = Integrals.adaptiveGaussianQuadrature(fx, Math.PI / 4.0, Math.PI / 2.0, 1e-12, 1e-12, 150);
		assertEquals(0.5832680936916835, qadrat, 1e-12);
	}
}



