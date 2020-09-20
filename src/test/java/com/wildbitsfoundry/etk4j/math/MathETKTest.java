package com.wildbitsfoundry.etk4j.math;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Test;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import static com.wildbitsfoundry.etk4j.math.MathETK.*;

public final class MathETKTest {

	@Test
	public void testQuadraticFormula() {
		Complex[] q = Formulas.quadraticFormula(2.0, 2.0, 1.0);
		assertArrayEquals(new Complex[] { new Complex(-0.5, 0.5), new Complex(-0.5, -0.5) }, q);
		
		q = Formulas.quadraticFormula(3.0, 12.0, 0.0);
		assertArrayEquals(new Complex[] { new Complex(0.0, 0.0), new Complex(-4.0, 0.0) }, q);
		
		q = Formulas.quadraticFormula(2.0, 0.0, 8.0);
		assertArrayEquals(new Complex[] { new Complex(0.0, 2.0), new Complex(0.0, -2.0) }, q);
		
		q = Formulas.quadraticFormula(2.0, 0.0, -8.0);
		assertArrayEquals(new Complex[] { new Complex(2.0, 0.0), new Complex(-2.0, 0.0) }, q);
		
		q = Formulas.quadraticFormula(2.0, 0.0, 0.0);
		assertArrayEquals(new Complex[] { new Complex(0.0, 0.0), new Complex(0.0, 0.0) }, q);
	}
	
	@Test
	public void testHypot() {
		double h = hypot(2.5, 9.4);
		assertEquals(9.726767191621274, h, 1e-12);
	}
	
	@Test
	public void testAsinh() {
		double a = asinh(2.5);
		assertEquals(1.647231146371096, a, 1e-12);
		
		a = asinh(-2.5);
		assertEquals(-1.647231146371096, a, 1e-12);
	}
	
	@Test
	public void testAcosh() {
		double a = acosh(2.5);
		assertEquals(1.566799236972411, a, 1e-12);
	}
	
	@Test
	public void testAtanh() {
		Complex a = atanh(2.5);
		assertEquals(0.423648930193602 , a.real(), 1e-12);
		assertEquals(1.570796326794897, a.imag(), 1e-12);
		
		a = atanh(-2.5);
		assertEquals(-0.423648930193602 , a.real(), 1e-12);
		assertEquals(-1.570796326794897, a.imag(), 1e-12);
		
		a = atanh(0.5);
		assertEquals(0.549306144334055 , a.real(), 1e-12);
		assertEquals(0.0, a.imag(), 1e-12);
		
		a = atanh(-0.5);
		assertEquals(-0.549306144334055 , a.real(), 1e-12);
		assertEquals(0.0, a.imag(), 1e-12);
	}
	
	@Test
	public void tetFix() {
		double f = 1.5;
		assertEquals(1.0 , fix(f), 1e-12);
		assertEquals(-1.0 , fix(-f), 1e-12);
	}
	
	@Test
	public void tetRem() {
		assertEquals(3.0 , rem(3, 4), 1e-12);
		assertEquals(3.0 , rem(3,-4), 1e-12);
		assertEquals(-3.0 , rem(-3,4), 1e-12);
		assertEquals(-3.0, rem(-3,-4), 1e-12);
		assertEquals(1.0 , rem(5,4), 1e-12);
		assertEquals(1.0 , rem(5,-4), 1e-12);
		assertEquals(-1.0 , rem(-5,4), 1e-12);
		assertEquals(-1.0 , rem(-5,-4), 1e-12);
		assertEquals(Double.NaN , rem(-5, 0), 1e-12);
	}

}
