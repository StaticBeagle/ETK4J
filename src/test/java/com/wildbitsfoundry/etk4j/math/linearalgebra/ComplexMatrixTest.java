package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;

public class ComplexMatrixTest {

    @Test
    public void testInverse() {
        Matrix A = new Matrix(new double[][]{{-2, -1}, {1, 0}});

        double w = 100;
        Complex jw = Complex.fromImaginary(w);
        ComplexMatrix cm = Matrix.identity(A.getRowCount()).multiply(jw).subtract(A);
        ComplexMatrix inv = cm.inv();
        double[] real = {1.999600059992001E-4, 9.99700049993001E-5, -9.99700049993001E-5, 1.9996000599920008E-8};
        double[] imag = {-0.00999700049993001, 1.9996000599920007E-6, -1.9996000599920007E-6, -0.010000999700049992};

        Complex[] expected = ComplexArrays.zip(real, imag);
        assertArrayEquals(expected, inv.getArray());

    }
}
