package com.wildbitsfoundry.etk4j.math.polynomials;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

public class PolynomialTest {
	
	@Test
	public void testConstructors() {
		Polynomial p1 = new Polynomial(0, 0, 0, 0, 0, 1, 2, 1);
		assertArrayEquals(new double[] { 1, 2, 1 }, p1.getCoefficients(), 1e-12);
		
		p1 = new Polynomial(0, 0, 0, 0, 0, 0);
		assertArrayEquals(new double[] { 0 }, p1.getCoefficients(), 1e-12);
		
		Polynomial p2 = new Polynomial(1, 2, 1);
		p2.getRoots();
		p1 = new Polynomial(p2);
		assertArrayEquals(new double[] { 1, 2, 1 }, p1.getCoefficients(), 1e-12);
	}

	@Test
	public void testGetDegree() {
		Polynomial poly = new Polynomial(5, 25, 50);
		int degree = poly.degree();
		assertEquals(2, degree, 1e-12);
	}

	@Test
	public void testNormalize() {
		Polynomial poly = new Polynomial(5, 25, 50);
		poly.normalize();
		assertArrayEquals(new double[] { 1, 5, 10 }, poly.getCoefficients(), 1e-12);
	}

	@Test
	public void testDenormalize() {
		Polynomial poly = new Polynomial(50, 25, 5);
		poly.denormalize();
		assertArrayEquals(new double[] { 10, 5, 1 }, poly.getCoefficients(), 1e-12);
		
		poly = new Polynomial(4, 0, 4, 0);
		poly.denormalize();
		assertArrayEquals(new double[] { 1, 0, 1, 0 }, poly.getCoefficients(), 1e-12);
		
		poly = new Polynomial(4, 0, 0, 0);
		poly.denormalize();
		assertArrayEquals(new double[] { 1, 0, 0, 0 }, poly.getCoefficients(), 1e-12);
	}

	@Test
	public void testMultiply() {
		Polynomial p1 = new Polynomial(1, 2, 1);
		Polynomial p2 = new Polynomial(-1, 4, 5, 6);

		Polynomial p3 = p1.multiply(p2);
		assertArrayEquals(new double[] { -1, 2, 12, 20, 17, 6 }, p3.getCoefficients(), 1e-12);

		p3 = p2.multiply(p1);
		assertArrayEquals(new double[] { -1, 2, 12, 20, 17, 6 }, p3.getCoefficients(), 1e-12);

		p3 = p1.multiply(-1.0, 4.0, 5.0, 6.0);
		assertArrayEquals(new double[] { -1, 2, 12, 20, 17, 6 }, p3.getCoefficients(), 1e-12);

		p3 = p1.multiply(5.0);
		assertArrayEquals(new double[] { 5, 10, 5 }, p3.getCoefficients(), 1e-12);
	}

	@Test
	public void testMultiplyEquals() {
		Polynomial p1 = new Polynomial(1, 2, 1);
		Polynomial p2 = new Polynomial(-1, 4, 5, 6);

		p1.multiplyEquals(p2);
		assertArrayEquals(new double[] { -1, 2, 12, 20, 17, 6 }, p1.getCoefficients(), 1e-12);

		p1.multiplyEquals(-9, 0, 4);
		assertArrayEquals(new double[] { 9, -18, -112, -172, -105, 26, 68, 24 }, p1.getCoefficients(), 1e-12);

		p1.multiplyEquals(-0.5);
		assertArrayEquals(new double[] { -4.5, 9, 56, 86, 52.5, -13, -34, -12 }, p1.getCoefficients(), 1e-12);
	}

	@Test
	public void testAdd() {
		Polynomial p1 = new Polynomial(1, 2, 1);
		Polynomial p2 = new Polynomial(-1, 4, -5, 6);

		Polynomial p3 = p1.add(p2);
		assertArrayEquals(new double[] { -1, 5, -3, 7 }, p3.getCoefficients(), 1e-12);

		p3 = p2.add(p1);
		assertArrayEquals(new double[] { -1, 5, -3, 7 }, p3.getCoefficients(), 1e-12);
	}

	@Test
	public void testAddEquals() {
		Polynomial p1 = new Polynomial(1, 2, 1);
		Polynomial p2 = new Polynomial(-1, 4, -5, 6);

		p1.addEquals(p2);
		assertArrayEquals(new double[] { -1, 5, -3, 7 }, p1.getCoefficients(), 1e-12);
	}

	@Test
	public void testSubtract() {
		Polynomial p1 = new Polynomial(1, 2, 1);
		Polynomial p2 = new Polynomial(-1, 4, -5, 6);

		Polynomial p3 = p1.subtract(p2);
		assertArrayEquals(new double[] { 1, -3, 7, -5 }, p3.getCoefficients(), 1e-12);

		p3 = p2.subtract(p1);
		assertArrayEquals(new double[] { -1, 3, -7, 5 }, p3.getCoefficients(), 1e-12);
	}

	@Test
	public void testSubtractEquals() {
		Polynomial p1 = new Polynomial(1, 2, 1);
		Polynomial p2 = new Polynomial(-1, 4, -5, 6);

		p1.subtractEquals(p2);
		assertArrayEquals(new double[] { 1, -3, 7, -5 }, p1.getCoefficients(), 1e-12);
	}

	@Test
	public void testPolyfit() {
		double[] x = { 1, 2, 3, 4 };
		double[] y = { 1, 10, 12, 15 };

		// Over determined
		Polynomial poly = Polynomial.polyFit(x, y, 2);
		assertArrayEquals(new double[] { -1.5, 11.9, -9.0 }, poly.getCoefficients(), 1e-12);

		// Unique solution
		Polynomial poly2 = Polynomial.polyFit(x, y, 3);
		double[] coefsUnique = new double[] { 1.333333333333, -11.5, 34.1666666666664, -23.0 };
		assertArrayEquals(coefsUnique, poly2.getCoefficients(), 1e-12);

		// Under determined
		Polynomial poly3 = Polynomial.polyFit(x, y, 4);
		double[] coefsUnder = new double[] { 0.6079433590792012, -4.746100257458709, 9.778017567772217,
				3.769498712706295, -8.40935938209905 };
		assertArrayEquals(coefsUnder, poly3.getCoefficients(), 1e-12);
	}

	@Test
	public void testDerivative() {
		Polynomial p1 = new Polynomial(1, 2, 1);
		Polynomial p2 = p1.derivative();

		assertArrayEquals(new double[] { 2, 2 }, p2.getCoefficients(), 1e-12);
	}
	
	@Test
	public void testDifferentiate() {
		Polynomial p1 = new Polynomial(1, 2, 1);
		double d = p1.differentiate(2);

		assertEquals(6, d, 1e-12);
	}
	
	@Test
	public void testIntegral() {
		Polynomial p1 = new Polynomial(3, 2, 1);
		Polynomial p2 = p1.integral(1.0);
		
		assertArrayEquals(new double[] { 1, 1, 1, 1 }, p2.getCoefficients(), 1e-12);
		
		p2 = p1.integral();
		
		assertArrayEquals(new double[] { 1, 1, 1, 0 }, p2.getCoefficients(), 1e-12);
	}

	@Test
	public void testIntegrate() {
		Polynomial p1 = new Polynomial(3, 2, 1);
		double d = p1.integrate(2, 4);
		
		assertEquals(70, d, 1e-12);
	}
	
	@Test
	public void testGetCoefficients() {
		Polynomial p1 = new Polynomial(1, 2, 1);

		assertArrayEquals(new double[] { 1, 2, 1 }, p1.getCoefficients(), 1e-12);
	}
	
	@Test
	public void testGetCoefficientAt() {
		Polynomial p1 = new Polynomial(1, 2, 1);

		assertEquals( 1, p1.getCoefficientAt(0), 1e-12);
	}

	@Test
	public void testEvaluateAt() {
		Polynomial p1 = new Polynomial(1, 2, 1);
		assertEquals(9.0, p1.evaluateAt(2), 1e-12);

		assertEquals(Complex.fromReal(9.0), p1.evaluateAt(2, 0));

		assertEquals(new Complex(5.0, 12.0), p1.evaluateAt(2, 2));

		assertEquals(new Complex(5.0, 12.0), p1.evaluateAt(new Complex(2, 2)));
	}
	
	@Test
	public void testReverseInPlace() {
		Polynomial p1 = new Polynomial(3, 2, 1);
		p1.reverseInPlace();
		assertArrayEquals(new double[] { 1, 2, 3 }, p1.getCoefficients(), 1e-12);
	}

	@Test
	public void testGetRoots() {
		double[] re = new double[] { -0.5, -0.5 };
		double[] im = new double[] { +0.5, -0.5 };
		Complex[] roots = ComplexArrays.zip(re, im);
		Polynomial p1 = new Polynomial(roots);

		double[] real = ComplexArrays.real(p1.getRoots());
		double[] imag = ComplexArrays.imag(p1.getRoots());

		assertArrayEquals(re, real, 1e-12);
		assertArrayEquals(im, imag, 1e-12);

		re = new double[] { -1, -2 };
		im = new double[] { +0, +0 };
		roots = ComplexArrays.zip(re, im);
		p1 = new Polynomial(roots);

		real = ComplexArrays.real(p1.getRoots());
		imag = ComplexArrays.imag(p1.getRoots());

		assertArrayEquals(re, real, 1e-12);
		assertArrayEquals(im, imag, 1e-12);

		re = new double[] { -1, -2 };
		im = new double[] { +0, +0 };
		roots = ComplexArrays.zip(re, im);
		p1 = new Polynomial(roots);

		real = ComplexArrays.real(p1.getRoots());
		imag = ComplexArrays.imag(p1.getRoots());

		assertArrayEquals(re, real, 1e-12);
		assertArrayEquals(im, imag, 1e-12);

		roots = new Polynomial(8, 6, 7, 5, 3, 0, 9).getRoots();
		re = new double[] { -0.9464031534859854, -0.9464031534859854, -0.061637902041647735, -0.061637902041647735,
				0.6330410555276348, 0.6330410555276348 };
		im = new double[] { 0.5578426205839989, -0.5578426205839989, 1.106706339734784, -1.106706339734784,
				0.5983158247369592, -0.5983158247369592 };

		real = ComplexArrays.real(roots);
		imag = ComplexArrays.imag(roots);

		assertArrayEquals(re, real, 1e-12);
		assertArrayEquals(im, imag, 1e-12);

		roots = new Polynomial(1, 0, 0, 0, 0, 0, 0).getRoots();
		re = new double[] { 0, 0, 0, 0, 0, 0 };
		im = new double[] { 0, 0, 0, 0, 0, 0 };

		real = ComplexArrays.real(roots);
		imag = ComplexArrays.imag(roots);

		assertArrayEquals(re, real, 1e-12);
		assertArrayEquals(im, imag, 1e-12);

		roots = new Polynomial(1, 1, 1, 1, 1, 1, 1).getRoots();
		re = new double[] { 0.623489801858734, 0.623489801858734, -0.9009688679024188, -0.9009688679024188,
				-0.22252093395631478, -0.22252093395631478 };
		im = new double[] { 0.7818314824680292, -0.7818314824680292, 0.43388373911755773, -0.43388373911755773,
				0.974927912181824, -0.974927912181824 };

		real = ComplexArrays.real(roots);
		imag = ComplexArrays.imag(roots);

		assertArrayEquals(re, real, 1e-12);
		assertArrayEquals(im, imag, 1e-12);
		
		roots = new Polynomial(1, 0, 0, 0.0, 1).getRoots();
		re = new double[] { -0.7071067811865477, -0.7071067811865477, 0.707106781186547, 0.707106781186547 };
		im = new double[] { 0.7071067811865478, -0.7071067811865478, 0.707106781186547, -0.707106781186547 };

		real = ComplexArrays.real(roots);
		imag = ComplexArrays.imag(roots);

		assertArrayEquals(re, real, 1e-12);
		assertArrayEquals(im, imag, 1e-12);
	
		roots = new Polynomial(2, 2, 1).getRoots();
		re = new double[] { -0.5, -0.5 };
		im = new double[] { 0.5, -0.5 };

		real = ComplexArrays.real(roots);
		imag = ComplexArrays.imag(roots);

		assertArrayEquals(re, real, 1e-12);
		assertArrayEquals(im, imag, 1e-12);
		
		roots = new Polynomial(1, 2, 1).getRoots();
		re = new double[] { -1, -1 };
		im = new double[] { 0, 0 };

		real = ComplexArrays.real(roots);
		imag = ComplexArrays.imag(roots);

		assertArrayEquals(re, real, 1e-12);
		assertArrayEquals(im, imag, 1e-12);
		
		roots = new Polynomial(1, 1).getRoots();
		re = new double[] { -1 };
		im = new double[] { 0 };

		real = ComplexArrays.real(roots);
		imag = ComplexArrays.imag(roots);

		assertArrayEquals(re, real, 1e-12);
		assertArrayEquals(im, imag, 1e-12);
		
		roots = new Polynomial().getRoots();

		assertEquals(0, roots.length);		
	}
	
	@Test
	public void testSubstite() {
		Polynomial p1 = new Polynomial(1, 2, 1);
		Polynomial p2 = p1.substitute(2);
		
		assertArrayEquals(new double[] { 4, 4, 1 }, p2.getCoefficients(), 1e-12);
		
		p1 = new Polynomial(1, 2, 1);
		p2 = p1.substitute(0.5);
		
		assertArrayEquals(new double[] { 0.25, 1, 1 }, p2.getCoefficients(), 1e-12);
		
		p2 = p1.substitute(0);
		
		assertArrayEquals(new double[] { 1 }, p2.getCoefficients(), 1e-12);
		
		p1.substituteInPlace(2);
		
		assertArrayEquals(new double[] { 4, 4, 1 }, p1.getCoefficients(), 1e-12);
		
		p1 = new Polynomial(1, 2, 1);
		p2 = p1.substitute(p1);
		
		assertArrayEquals(new double[] { 1, 4, 8, 8, 4 }, p2.getCoefficients(), 1e-12);
		
		p1 = new Polynomial(1, 0, -1, 2, -1);
		p2 = p1.substitute(new Polynomial(1, 2, 1));
		
		assertArrayEquals(new double[] { 1, 8, 28, 56, 69, 52, 24, 8, 1 }, p2.getCoefficients(), 1e-12);
		
		RationalFunction rf = new Polynomial(1, 1, 2, 1).substitute(new Polynomial(1, 0), new Polynomial(1, 2, 1));
		
		p1 = rf.getNumerator();
		p2 = rf.getDenominator();
		
		assertArrayEquals(new double[] { 1, 8, 24, 35, 24, 8, 1 }, p1.getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 6, 15, 20, 15, 6, 1}, p2.getCoefficients(), 1e-12);
	}
	
	@Test
	public void testSubstituteInPlace() {
		Polynomial p1 = new Polynomial(1, 2, 1);
		p1.substituteInPlace(p1);
		
		assertArrayEquals(new double[] { 1, 4, 8, 8, 4 }, p1.getCoefficients(), 1e-12);
	}
	
	@Rule
	public final ExpectedException exception = ExpectedException.none();
	
	
	@Test
	public void testPow() {
		
		Polynomial p1 = new Polynomial(1, 2, 1);
		Polynomial p2 = p1.pow(0);
		
		assertArrayEquals(new double[] { 1 }, p2.getCoefficients(), 1e-12);
		
		p2 = p1.pow(3);
		
		assertArrayEquals(new double[] { 1, 6, 15, 20, 15, 6, 1 }, p2.getCoefficients(), 1e-12);
		
		exception.expect(IllegalArgumentException.class);
		p1.pow(-1);
	}
	
	@Test
	public void testToString() {
		Polynomial p1 = new Polynomial(1, 2, 1);
		
		assertEquals("x^2 + 2.000 * x + 1.000", p1.toString());
		
		p1 = new Polynomial(-1, 2, -1, 0, 1, 1);
		assertEquals("-x^5 + 2.000 * x^4 - x^3 + x + 1.000", p1.toString());
		
		p1 = new Polynomial(1);
		assertEquals("1.000", p1.toString());
		
		p1 = new Polynomial(1, 0);
		assertEquals("x", p1.toString());
		
		p1 = new Polynomial(0);
		assertEquals("0.000", p1.toString());
	}
	
	
}
