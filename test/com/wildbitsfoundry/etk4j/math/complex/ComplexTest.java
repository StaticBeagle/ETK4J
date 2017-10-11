package com.wildbitsfoundry.etk4j.math.complex;

import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public class ComplexTest {
	
	Complex a;
	Complex b;
	
	@Before
	public void before() {
		a = new Complex(-2, 3);
		b = new Complex(5, -6);
	}
	
	@Test
	public void testAbs() {
		assertEquals(3.605551275463989, a.abs(), 1e-12);
	}
	
	@Test
	public void testArg() {
		assertEquals( 2.158798930342464, a.arg(), 1e-12);
	}
	
	@Test
	public void testNorm() {
		assertEquals(13.0, a.norm(),  1e-12);
	}
	
	@Test
	public void testInvert() {
		Complex c = a.invert();
		assertEquals(-0.153846153846154, c.real(), 1e-12);
		assertEquals(-0.230769230769231, c.imag(), 1e-12);
	}
	
	@Test
	public void testAdd() {
		Complex c = a.add(b);
		assertEquals(3.0, c.real(), 1e-12);
		assertEquals(-3.0, c.imag(), 1e-12);
		
		c = a.add(-1);
		assertEquals(-3.0, c.real(), 1e-12);
		assertEquals(3.0, c.imag(), 1e-12);
		
		c.addEquals(b);
		assertEquals(2.0, c.real(), 1e-12);
		assertEquals(-3.0, c.imag(), 1e-12);
		
		c.addEquals(1);
		assertEquals(3.0, c.real(), 1e-12);
		assertEquals(-3.0, c.imag(), 1e-12);
	}
	
	@Test
	public void testSubtract() {
		Complex c = a.subtract(b);
		assertEquals(-7.0, c.real(), 1e-12);
		assertEquals(9.0, c.imag(), 1e-12);
		
		c = a.subtract(-1);
		assertEquals(-1.0, c.real(), 1e-12);
		assertEquals(3.0, c.imag(), 1e-12);
		
		c.subtractEquals(b);
		assertEquals(-6.0, c.real(), 1e-12);
		assertEquals(9.0, c.imag(), 1e-12);
		
		c.addEquals(1);
		assertEquals(-5.0, c.real(), 1e-12);
		assertEquals(9.0, c.imag(), 1e-12);
	}
	
	@Test
	public void testMultiply() {
		Complex c = a.multiply(b);
		assertEquals(8.0, c.real(), 1e-12);
		assertEquals(27.0, c.imag(), 1e-12);
		
		c = a.multiply(2);
		assertEquals(-4.0, c.real(), 1e-12);
		assertEquals(6.0, c.imag(), 1e-12);
		
		c.multiplyEquals(b);
		assertEquals(16.0, c.real(), 1e-12);
		assertEquals(54.0, c.imag(), 1e-12);
		
		c.multiplyEquals(2);
		assertEquals(32.0, c.real(), 1e-12);
		assertEquals(108.0, c.imag(), 1e-12);
	}
	
	@Test
	public void testDivide() {
		Complex c = a.divide(b);
		assertEquals(-0.459016393442623, c.real(), 1e-12);
		assertEquals(0.049180327868852, c.imag(), 1e-12);
		
		c = a.divide(2);
		assertEquals(-1.0, c.real(), 1e-12);
		assertEquals(1.5, c.imag(), 1e-12);
		
		c.divideEquals(b);
		assertEquals(-0.229508196721311, c.real(), 1e-12);
		assertEquals(0.024590163934426, c.imag(), 1e-12);

		c.divideEquals(2);
		assertEquals(-0.114754098360656, c.real(), 1e-12);
		assertEquals(0.012295081967213, c.imag(), 1e-12);
	}
	
	@Test
	public void testSqrt() {
		Complex c = a.sqrt();
		assertEquals(0.8959774761298381, c.real(), 1e-12);
		assertEquals(1.6741492280355401, c.imag(), 1e-12);
	}
	
	@Test
	public void testPow() {
		Complex c = a.pow(b);
		assertEquals(-2.5691757471583375E8, c.real(), 1e-12);
		assertEquals(1.091169633199824E7, c.imag(), 1e-12);
				
		c = a.pow(-2.0);
		assertEquals(-0.029585798816568018, c.real(), 1e-12);
		assertEquals(0.07100591715976333, c.imag(), 1e-12);
		
		c = a.pow(-1.0);
		assertEquals(-0.15384615384615388, c.real(), 1e-12);
		assertEquals(-0.23076923076923075, c.imag(), 1e-12);
		
		c = a.pow(-0.5);
		assertEquals(0.2484994409113033, c.real(), 1e-12);
		assertEquals(-0.46432545265081493, c.imag(), 1e-12);

		c = a.pow(0.0);
		assertEquals(1.0, c.real(), 1e-12);
		assertEquals(0.0, c.imag(), 1e-12);
		
		c = a.pow(0.5);
		assertEquals(0.895977476129838, c.real(), 1e-12);
		assertEquals(1.6741492280355401, c.imag(), 1e-12);
		
		c = a.pow(1.0);
		assertEquals(-2.0, c.real(), 1e-12);
		assertEquals(3.0, c.imag(), 1e-12);
		
		c = a.pow(2.0);
		assertEquals(-5.0, c.real(), 1e-12);
		assertEquals(-12.0, c.imag(), 1e-12);
	}
	
	@Test
	public void testLog() {
		Complex c = a.log();
		assertEquals(1.2824746787307684, c.real(), 1e-12);
		assertEquals(2.1587989303424644, c.imag(), 1e-12);
	}
	
	@Test
	public void testExp() {
		Complex c = a.exp();
		assertEquals(-0.13398091492954262, c.real(), 1e-12);
		assertEquals(0.019098516261135196, c.imag(), 1e-12);
	}
	
	@Test
	public void testSin() {
		Complex c = a.sin();
		assertEquals(-9.15449914691143, c.real(), 1e-12);
		assertEquals(-4.168906959966565, c.imag(), 1e-12);
	}
}
