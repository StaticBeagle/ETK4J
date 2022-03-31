package com.wildbitsfoundry.etk4j.util;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;

public class ComplexArraysTest {

    @Test
    public void testAddition() {
        Complex[] a = {new Complex(1,1), new Complex(2, 2), new Complex(3, 3)};
        Complex b = new Complex(4, 4);
        Complex[] expected = {new Complex(5, 5), new Complex(6, 6), new Complex(7 ,7)};
        assertArrayEquals(expected, ComplexArrays.addElementWise(a, b));
    }

    @Test
    public void testAdditionArrays() {
        Complex[] a = {new Complex(1,1), new Complex(2, 2), new Complex(3, 3)};
        Complex[] b = a;
        Complex[] expected = {new Complex(2, 2), new Complex(4, 4), new Complex(6 ,6)};
        assertArrayEquals(expected, ComplexArrays.addElementWise(a, b));
    }

    @Test
    public void multiplyElementWise() {
        Complex[] a = {new Complex(1,1), new Complex(2, 2), new Complex(3, 3)};
        double b = 4;
        Complex[] expected = {new Complex(4, 4), new Complex(8, 8), new Complex(12 ,12)};
        assertArrayEquals(expected, ComplexArrays.multiplyElementWise(a, b));
    }

    @Test
    public void testSubtraction() {
        Complex[] a = {new Complex(1,1), new Complex(2, 2), new Complex(3, 3)};
        Complex b = new Complex(4, 4);
        Complex[] expected = {new Complex(-3, -3), new Complex(-2, -2), new Complex(-1 ,-1)};
        assertArrayEquals(expected, ComplexArrays.subtractElementWise(a, b));
    }

    @Test
    public void testSubtractionRealValue() {
        Complex[] a = {new Complex(1,1), new Complex(2, 2), new Complex(3, 3)};
        double b = 4;
        Complex[] expected = {new Complex(-3, 1), new Complex(-2, 2), new Complex(-1 ,3)};
        assertArrayEquals(expected, ComplexArrays.subtractElementWise(a, b));
    }
}
