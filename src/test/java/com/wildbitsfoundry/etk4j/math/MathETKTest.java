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
		assertArrayEquals(new Complex[] { new Complex(), Complex.fromReal(-4.0) }, q);
		
		q = Formulas.quadraticFormula(2.0, 0.0, 8.0);
		assertArrayEquals(new Complex[] { Complex.fromImaginary(2.0), Complex.fromImaginary(-2.0) }, q);
		
		q = Formulas.quadraticFormula(2.0, 0.0, -8.0);
		assertArrayEquals(new Complex[] { Complex.fromReal(2.0), Complex.fromReal(-2.0) }, q);
		
		q = Formulas.quadraticFormula(2.0, 0.0, 0.0);
		assertArrayEquals(new Complex[] { new Complex(), new Complex() }, q);

		q = Formulas.quadraticFormula(1.0, 5.0, 6.0);
		assertArrayEquals(new Complex[] { Complex.fromReal(-3.0), Complex.fromReal(-2.0) }, q);
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

	@Test
	public void testFrexp() {
		double value = -5.35;
		MathETK.FRexpResult frexp = frexp(value);
		assertEquals(frexp.exponent, 3);
		assertEquals(frexp.mantissa, -0.66875, 1e-12);
		assertEquals(-5.35, frexp.mantissa * (1 << frexp.exponent), 1e-12);

		value = 5.35;
		frexp = frexp(value);
		assertEquals(frexp.exponent, 3);
		assertEquals(frexp.mantissa, 0.66875, 1e-12);
		assertEquals(5.35, frexp.mantissa * (1 << frexp.exponent), 1e-12);

		value = 8.0;
		frexp = frexp(value);
		assertEquals(frexp.exponent, 3);
		assertEquals(frexp.mantissa, 1.0, 1e-12);
		assertEquals(8.0, frexp.mantissa * (1 << frexp.exponent), 1e-12);

		value = 0.0;
		frexp = frexp(value);
		assertEquals(frexp.exponent, 0);
		assertEquals(frexp.mantissa, 0.0, 1e-12);
		assertEquals(0.0, frexp.mantissa * (1 << frexp.exponent), 1e-12);

		value = -0.0;
		frexp = frexp(value);
		assertEquals(frexp.exponent, 0);
		assertEquals(frexp.mantissa, 0.0, 1e-12);
		assertEquals(0.0, frexp.mantissa * (1 << frexp.exponent), 1e-12);

		value = Double.NaN;
		frexp = frexp(value);
		assertEquals(frexp.exponent, -1);
		assertEquals(frexp.mantissa, Double.NaN, 1e-12);
		assertEquals(Double.NaN, frexp.mantissa * (1 << frexp.exponent), 1e-12);

		value = Double.NEGATIVE_INFINITY;
		frexp = frexp(value);
		assertEquals(frexp.exponent, -1);
		assertEquals(frexp.mantissa, Double.NEGATIVE_INFINITY, 1e-12);
		assertEquals(Double.POSITIVE_INFINITY, frexp.mantissa * (1 << frexp.exponent), 1e-12);

		value = Double.POSITIVE_INFINITY;
		frexp = frexp(value);
		assertEquals(frexp.exponent, -1);
		assertEquals(frexp.mantissa, Double.POSITIVE_INFINITY, 1e-12);
		assertEquals(Double.NEGATIVE_INFINITY, frexp.mantissa * (1 << frexp.exponent), 1e-12);
	}

}
