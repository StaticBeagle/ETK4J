package com.wildbitsfoundry.etk4j.math.polynomial;

import org.junit.Test;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

public class RationalFunctionTest {
	
	@Test
	public void testConstructors() {
		Complex[] zeros = { Complex.fromReal(-1.0) };
		Complex[] poles = { Complex.fromReal(-1.0), Complex.fromReal(-1.0) };
		
		RationalFunction rf = new RationalFunction(zeros, poles);
		
		assertArrayEquals(new double[] { 1, 1 }, rf.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 2, 1 }, rf.getDenominator().getCoefficients(), 1e-12);
		
		rf = new RationalFunction(new Polynomial(1), new Polynomial(poles));
		assertArrayEquals(new double[] { 1 }, rf.getNumerator().getCoefficients(), 1e-12);
		
		RationalFunction rf2 = new RationalFunction(rf);
		
		assertArrayEquals(new double[] { 1 }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 2, 1 }, rf2.getDenominator().getCoefficients(), 1e-12);
		
		rf = new RationalFunction(new Polynomial(1, 1), new Polynomial(1, 2, 1));
		
		assertArrayEquals(new double[] { 1, 1 }, rf.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 2, 1 }, rf.getDenominator().getCoefficients(), 1e-12);
		
		rf = new RationalFunction(new Polynomial(1), new Polynomial(1, 2, 1));
		
		assertArrayEquals(new double[] { 1 }, rf.getNumerator().getCoefficients(), 1e-12);
		
		rf = new RationalFunction(new double[] {1, 1}, new double[] {1, 2, 1});
		
		assertArrayEquals(new double[] { 1 }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 2, 1 }, rf2.getDenominator().getCoefficients(), 1e-12);
	}
	
	@Test
	public void testGetZeros() {
		RationalFunction rf = new RationalFunction(new double[] {1, 1}, new double[] {1, 2, 1});
		Complex[] zeros = rf.getZeros();
		
		double[] real = ComplexArrays.real(zeros);
		double[] imag = ComplexArrays.imag(zeros);
		assertArrayEquals(new double[] { -1 }, real, 1e-12);
		assertArrayEquals(new double[] { 0 }, imag, 1e-12);
	}
	
	@Test
	public void testGetPoles() {
		RationalFunction rf = new RationalFunction(new double[] {1, 1}, new double[] {1, 2, 1});
		Complex[] zeros = rf.getPoles();
		
		double[] real = ComplexArrays.real(zeros);
		double[] imag = ComplexArrays.imag(zeros);
		assertArrayEquals(new double[] { -1, -1 }, real, 1e-12);
		assertArrayEquals(new double[] { 0, 0 }, imag, 1e-12);
	}
	
	@Test
	public void testAdd() {
		RationalFunction rf = new RationalFunction(new double[] {1, 1}, new double[] {1, 2, 1});
		RationalFunction rf2 = rf.add(rf);
		
		assertArrayEquals(new double[] { 2, 6, 6, 2 }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 4, 6, 4, 1 }, rf2.getDenominator().getCoefficients(), 1e-12);
		
		rf2 = rf.add(2);
		
		assertArrayEquals(new double[] { 2, 5, 3 }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 2, 1 }, rf2.getDenominator().getCoefficients(), 1e-12);
	}
	
	@Test
	public void testSubtract() {
		RationalFunction rf = new RationalFunction(new double[] {1, 1}, new double[] {1, 2, 1});
		RationalFunction rf2 = rf.subtract(rf);
		
		assertArrayEquals(new double[] { 0 }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 4, 6, 4, 1 }, rf2.getDenominator().getCoefficients(), 1e-12);
		
		rf2 = rf.subtract(rf.multiply(2));
		
		assertArrayEquals(new double[] { -1, -3, -3, -1 }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 4, 6, 4, 1 }, rf2.getDenominator().getCoefficients(), 1e-12);
	}
	
	@Test
	public void testMultiply() {
		RationalFunction rf = new RationalFunction(new double[] {1, 1}, new double[] {1, 2, 1});
		RationalFunction rf2 = rf.multiply(2);
		
		assertArrayEquals(new double[] { 2, 2 }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 2, 1 }, rf2.getDenominator().getCoefficients(), 1e-12);
		
		rf2 = rf.multiply(new Polynomial(1, 1));
		
		assertArrayEquals(new double[] { 1, 2, 1  }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 2, 1  }, rf2.getDenominator().getCoefficients(), 1e-12);
		
		rf2 = rf.multiply(rf);
		
		assertArrayEquals(new double[] { 1, 2, 1  }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 4, 6, 4, 1 }, rf2.getDenominator().getCoefficients(), 1e-12);
	}
	
	@Test
	public void testSubstitute() {
		RationalFunction rf = new RationalFunction(new double[] {1, 1}, new double[] {1, 2, 1});
		RationalFunction rf2 = rf.substitute(2);
		
		assertArrayEquals(new double[] { 2, 1 }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 4, 4, 1 }, rf2.getDenominator().getCoefficients(), 1e-12);
		
		rf2.substituteInPlace(0);
		
		assertArrayEquals(new double[] { 1 }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1 }, rf2.getDenominator().getCoefficients(), 1e-12);
		
		rf2 = rf.substitute(new Polynomial(1, 1));
		
		assertArrayEquals(new double[] { 1, 2 }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 4, 4 }, rf2.getDenominator().getCoefficients(), 1e-12);
		
		rf2.substituteInPlace(new Polynomial(1, 1));
		
		assertArrayEquals(new double[] { 1, 3 }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 6, 9 }, rf2.getDenominator().getCoefficients(), 1e-12);
		
		rf2 = rf.substitute(new RationalFunction(new Polynomial(1, 0), new Polynomial(1, 1)));
		
		assertArrayEquals(new double[] { 2, 3, 1 }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 4, 4, 1 }, rf2.getDenominator().getCoefficients(), 1e-12);
		
		rf = new RationalFunction(new double[] {1, 2, 1}, new double[] {1, 1});
		rf2 = rf.substitute(new RationalFunction(new Polynomial(1, 0), new Polynomial(1, 1)));
		
		assertArrayEquals(new double[] { 4, 4, 1 }, rf2.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 2, 3, 1 }, rf2.getDenominator().getCoefficients(), 1e-12);
		
		rf.substituteInPlace(new RationalFunction(new Polynomial(1, 0), new Polynomial(1, 1)));
		
		assertArrayEquals(new double[] { 4, 4, 1 }, rf.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 2, 3, 1 }, rf.getDenominator().getCoefficients(), 1e-12);
	}
	
	@Test
	public void testNormalize() {
		RationalFunction rf = new RationalFunction(new double[] {2, 2}, new double[] {2, 4, 2});
		double norm = rf.normalize();
		
		assertEquals(1.0, norm, 1e-12);
		assertArrayEquals(new double[] { 1, 1 }, rf.getNumerator().getCoefficients(), 1e-12);
		assertArrayEquals(new double[] { 1, 2, 1 }, rf.getDenominator().getCoefficients(), 1e-12);
	}
	
	@Test
	public void testEvaluate() {
		RationalFunction rf = new RationalFunction(new double[] {1, 1}, new double[] {1, 2, 1});
		double val = rf.evaluateAt(0.0);
		
		assertEquals(1.0, val, 1e-12);

		val = rf.evaluateAt(2.0);
		
		assertEquals(1.0 / 3.0, val, 1e-12);
		
		Complex cval = rf.evaluateAt(0, 0);
		
		assertEquals(1.0, cval.real(), 1e-12);
		assertEquals(0.0, cval.imag(), 1e-12);
		
		cval = rf.evaluateAt(1, 1);
		
		assertEquals(0.4, cval.real(), 1e-12);
		assertEquals(-0.2, cval.imag(), 1e-12);
	}

}
