package com.wildbitsfoundry.etk4j.math.complex;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import org.junit.Test;

public class ComplexTest {


	@Test
	public void testHashCode() {
		Complex a = new Complex(-2.0, 3.0);
		new Complex(5.0, -6.0);
		int hash = a.hashCode();
		
		assertEquals(-2131229759, hash);
	}
	
	@Test
	public void testEquals() {
		Complex a = new Complex(-2.0, 3.0);
		new Complex(5.0, -6.0);

		Complex c = new Complex(a.real(), a.imag());
		assertEquals(a, c);
		
		Complex d = new Complex(a);
		assertEquals(a, d);
		
		c = a;
		assertEquals(a, c);
		
		assertNotEquals(a, c.toString());

		c = null;
		assertNotEquals(a, c);
		
		c = a.add(1);
		assertNotEquals(a, c);
		
		c = a.add(Complex.fromImaginary(1.0));
		assertNotEquals(a, c);
	}
	
	@Test
	public void testCompareTo() {
		Complex a = new Complex(-2.0, 3.0);
		new Complex(5.0, -6.0);
		Complex c = new Complex(a.real(), a.imag());
		
		assertEquals(0.0, a.compareTo(c), 1e-12);
		assertEquals(0.0, a.compareTo(a), 1e-12);
		assertEquals(1.0, a.compareTo(c.subtract(1.0)), 1e-12);
		assertEquals(-1.0, a.compareTo(c.add(1.0)), 1e-12);
	}
	
	@Test
	public void testFromPolar() {
		Complex c = Complex.fromPolar(0.5, 1.5);
		assertEquals(0.035368600833851, c.real(), 1e-12);
		assertEquals(0.498747493302027, c.imag(), 1e-12);
	}
	
	@Test
	public void testConj() {
		Complex a = new Complex(-2.0, 3.0);
		Complex c = a.conj();
		assertEquals(a.real(), c.real(), 1e-12);
		assertEquals(-a.imag(), c.imag(), 1e-12);
	}
	
	@Test
	public void testAbs() {
		Complex a = new Complex(-2.0, 3.0);
		assertEquals(3.605551275463989, a.abs(), 1e-12);
	}
	
	@Test
	public void testArg() {
		Complex a = new Complex(-2.0, 3.0);
		assertEquals( 2.158798930342464, a.arg(), 1e-12);
	}
	
	@Test
	public void testNorm() {
		Complex a = new Complex(-2.0, 3.0);
		assertEquals(13.0, a.norm(),  1e-12);
	}
	
	@Test
	public void testInvert() {
		Complex a = new Complex(-2.0, 3.0);
		Complex c = a.invert();
		assertEquals(-0.153846153846154, c.real(), 1e-12);
		assertEquals(-0.230769230769231, c.imag(), 1e-12);
		
		c = a;
		c.invertEquals();
		assertEquals(-0.153846153846154, c.real(), 1e-12);
		assertEquals(-0.230769230769231, c.imag(), 1e-12);
	}
	
	@Test
	public void testAdd() {
		Complex a = new Complex(-2.0, 3.0);
		Complex b = new Complex(5.0, -6.0);

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
		Complex a = new Complex(-2.0, 3.0);
		Complex b = new Complex(5.0, -6.0);

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
		
		c.subtractEquals(1);
		assertEquals(-6.0, c.real(), 1e-12);
		assertEquals(9.0, c.imag(), 1e-12);
	}
	
	@Test
	public void testMultiply() {
		Complex a = new Complex(-2.0, 3.0);
		Complex b = new Complex(5.0, -6.0);

		Complex c = a.multiply(b);
		assertEquals(8.0, c.real(), 1e-12);
		assertEquals(27.0, c.imag(), 1e-12);
		
		c = a.multiply(b.real(), b.imag());
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
		Complex a = new Complex(-2.0, 3.0);
		Complex b = new Complex(5.0, -6.0);

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
		Complex a = new Complex(-2.0, 3.0);
		new Complex(5.0, -6.0);

		Complex c = a.sqrt();
		assertEquals(0.8959774761298381, c.real(), 1e-12);
		assertEquals(1.6741492280355401, c.imag(), 1e-12);
		
		c = new Complex().sqrt();
		assertEquals(0.0, c.real(), 1e-12);
		assertEquals(0.0, c.imag(), 1e-12);
		
		c = a.uminus().sqrt();
		assertEquals(1.6741492280355401, c.real(), 1e-12);
		assertEquals(-0.8959774761298381, c.imag(), 1e-12);		
	}
	
	@Test
	public void testPow() {
		Complex a = new Complex(-2.0, 3.0);
		Complex b = new Complex(5.0, -6.0);

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
	public void testPow2() {
		Complex a = new Complex(-2.0, 3.0);
		new Complex(5.0, -6.0);

		Complex c = a.pow2();
		assertEquals(-5.0, c.real(), 1e-12);
		assertEquals(-12.0, c.imag(), 1e-12);
	}
	@Test
	public void testLog() {
		Complex a = new Complex(-2.0, 3.0);

		Complex c = a.log();
		assertEquals(1.2824746787307684, c.real(), 1e-12);
		assertEquals(2.1587989303424644, c.imag(), 1e-12);
	}
	
	@Test
	public void testExp() {
		Complex a = new Complex(-2.0, 3.0);

		Complex c = a.exp();
		assertEquals(-0.13398091492954262, c.real(), 1e-12);
		assertEquals(0.019098516261135196, c.imag(), 1e-12);
	}
	
	@Test
	public void testSin() {
		Complex a = new Complex(-2.0, 3.0);

		Complex c = a.sin();
		assertEquals(-9.15449914691143, c.real(), 1e-12);
		assertEquals(-4.168906959966565, c.imag(), 1e-12);
	}
	
	@Test
	public void testAsin() {
		Complex a = new Complex(-2.0, 3.0);

		Complex c = a.asin();
		assertEquals(-0.570652784321099, c.real(), 1e-12);
		assertEquals(1.983387029916536, c.imag(), 1e-12);
	}

	@Test
	public void testCos() {
		Complex a = new Complex(-2.0, 3.0);

		Complex c = a.cos();
		assertEquals(-4.189625690968807, c.real(), 1e-12);
		assertEquals(9.109227893755337, c.imag(), 1e-12);
	}
	
	@Test
	public void testAcos() {
		Complex a = new Complex(-2.0, 3.0);

		Complex c = a.acos();
		assertEquals(2.141449111115996, c.real(), 1e-12);
		assertEquals(-1.983387029916536, c.imag(), 1e-12);
	}

	@Test
	public void testTan() {
		Complex a = new Complex(-2.0, 3.0);

		Complex c = a.tan();
		assertEquals(0.003764025641504, c.real(), 1e-12);
		assertEquals(1.003238627353610, c.imag(), 1e-12);
		
		c = new Complex(1, 20.5).tan();
		assertEquals(0.0, c.real(), 1e-12);
		assertEquals(1.0, c.imag(), 1e-12);
		
		c = new Complex(1, -20.5).tan();
		assertEquals(0.0, c.real(), 1e-12);
		assertEquals(-1.00, c.imag(), 1e-12);
	}
	
	@Test
	public void testAtan() {
		Complex a = new Complex(-2.0, 3.0);

		Complex c = a.atan();
		assertEquals(-1.409921049596576, c.real(), 1e-12);
		assertEquals(0.229072682968539, c.imag(), 1e-12);
	}

	@Test
	public void testSinh() {
		Complex a = new Complex(-2.0, 3.0);

		Complex c = a.sinh();
		assertEquals(3.59056458998578, c.real(), 1e-12);
		assertEquals(0.5309210862485197, c.imag(), 1e-12);
	}

	@Test
	public void testCosh() {
		Complex a = new Complex(-2.0, 3.0);

		Complex c = a.cosh();
		assertEquals(-3.7245455049153224, c.real(), 1e-12);
		assertEquals(-0.5118225699873846, c.imag(), 1e-12);
	}
	
	@Test
	public void testTanh() {
		Complex a = new Complex(-2.0, 3.0);

		Complex c = a.tanh();
		assertEquals(-0.965385879022133, c.real(), 1e-12);
		assertEquals(-0.009884375038323, c.imag(), 1e-12);
		
		c = new Complex(1, 20.5).tanh();
		assertEquals(1.3070443537632264, c.real(), 1e-12);
		assertEquals(-0.05716427993519685, c.imag(), 1e-12);
		
		c = new Complex(1, -20.5).tanh();
		assertEquals(1.3070443537632264, c.real(), 1e-12);
		assertEquals(0.05716427993519685, c.imag(), 1e-12);
	}
	
	
	@Test
	public void testToString() {
		Complex a = new Complex(-2.0, 3.0);
		assertEquals("(-2.0000 + 3.0000j)", a.toString());
	}
}
