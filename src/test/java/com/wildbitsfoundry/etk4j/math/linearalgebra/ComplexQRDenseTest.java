package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;

public class ComplexQRDenseTest {

    @Test
    public void testGetQ() {
        // Define a 2x2 complex matrix
        Complex[][] A = {
                {new Complex(1, 1), new Complex(2, 2)},
                {new Complex(3, 3), new Complex(4, 4)}
        };

        // Perform Householder QR decomposition
        ComplexQRDecompositionDense complexQRDecompositionDense = new ComplexQRDecompositionDense(ComplexMatrixDense.from2DArray(A));
        Complex[] expected = {new Complex(-0.3162277660168378, 0), new Complex(0.9486832980505124, 0),
                new Complex(-0.9486832980505134, 0), new Complex(-0.3162277660168382, 0)};
        assertArrayEquals(expected, complexQRDecompositionDense.getQ().getArray());
    }

    @Test
    public void testGetR() {
        // Define a 2x2 complex matrix
        Complex[][] A = {
                {new Complex(1, 1), new Complex(2, 2)},
                {new Complex(3, 3), new Complex(4, 4)}
        };

        // Perform Householder QR decomposition
        ComplexQRDecompositionDense complexQRDecompositionDense = new ComplexQRDecompositionDense(ComplexMatrixDense.from2DArray(A));
        Complex[] expected = {new Complex(-3.1622776601683773, -3.1622776601683773), new Complex(-4.427188724235729, -4.427188724235729),
                new Complex(0, 0), new Complex(0.6324555320336729, 0.6324555320336729)};
        assertArrayEquals(expected, complexQRDecompositionDense.getR().getArray());
    }

    @Test
    public void testGetH() {
        // Define a 2x2 complex matrix
        Complex[][] A = {
                {new Complex(1, 1), new Complex(2, 2)},
                {new Complex(3, 3), new Complex(4, 4)}
        };

        // Perform Householder QR decomposition
        ComplexQRDecompositionDense complexQRDecompositionDense = new ComplexQRDecompositionDense(ComplexMatrixDense.from2DArray(A));
        Complex[] expected = {new Complex(0.8112421851755608, 0.8112421851755608), new Complex(0, 0),
                new Complex(0.5847102846637647, 0.5847102846637647), new Complex(-0.9999999999999998, -0.9999999999999998)};
        assertArrayEquals(expected, complexQRDecompositionDense.getH().getArray());
    }

    @Test
    public void testGetOriginalMatrixBack() {
        // Define a 2x2 complex matrix
        Complex[][] A = {
                {new Complex(1, 1), new Complex(2, 2)},
                {new Complex(3, 3), new Complex(4, 4)}
        };

        ComplexMatrixDense complexMatrixDense = ComplexMatrixDense.from2DArray(A);
        Complex[] expected = {new Complex(0.9999999999999989, 0.9999999999999989), new Complex(1.9999999999999951, 1.9999999999999951),
                new Complex(2.999999999999997, 2.999999999999997), new Complex(3.9999999999999973, 3.9999999999999973)};

        // Perform Householder QR decomposition
        ComplexQRDecompositionDense complexQRDecompositionDense = new ComplexQRDecompositionDense(complexMatrixDense);
        assertArrayEquals(expected, complexQRDecompositionDense.getQ().multiply(complexQRDecompositionDense.getR()).getArray());
    }

    @Test
    public void testGetQConjugateTranspose() {
        ComplexMatrixDense A = new ComplexMatrixDense(new Complex[][]{
                {new Complex(1, 1), new Complex(2, 0)},
                {new Complex(2, -1), new Complex(3, 0)}
        });

        // Perform Householder QR decomposition
        ComplexQRDecompositionDense complexQRDecompositionDense = new ComplexQRDecompositionDense(A);
        Complex[] expected = {new Complex(-0.5345224838248488, 0), new Complex(-0.26726124191242434, -0.8017837257372731),
                new Complex(0.26726124191242423, -0.8017837257372731), new Complex(-0.5345224838248488, 5.551115123125783E-17)};
        assertArrayEquals(expected, complexQRDecompositionDense.getQH().getArray());
    }

    @Test
    public void testQisOrthogonal() {
        // Define a 2x2 complex matrix
        Complex[][] A = {
                {new Complex(1, 1), new Complex(2, 2)},
                {new Complex(3, 3), new Complex(4, 4)}
        };

        Complex[] expected = {new Complex(0.9999999999999973, 0), new Complex(-1.1102230246251565E-16, 0),
                new Complex(-1.6653345369377348E-16, 0), new Complex(0.9999999999999996, 0)};

        // Perform Householder QR decomposition
        ComplexQRDecompositionDense complexQRDecompositionDense = new ComplexQRDecompositionDense(ComplexMatrixDense.from2DArray(A));
        assertArrayEquals(expected, complexQRDecompositionDense.getQ().multiply(complexQRDecompositionDense.getQH()).getArray());
    }

    @Test
    public void testSolve() {
        Complex[] x = {new Complex(1, 1), new Complex(2, 2)};
        Complex[] u = {new Complex(3, 3), new Complex(4, 4)};

        Complex[][] y = {x, u};
        Complex[] g = {new Complex(60, 70)};
        Complex[] f = {new Complex(80, 90)};

        Complex[] expected = {new Complex(-45.00000000000017, -5.000000000000018),
                new Complex(55.000000000000114, 5.000000000000011)};

        ComplexMatrixDense sol = new ComplexQRDecompositionDense(new ComplexMatrixDense(y))
                .solve(new ComplexMatrixDense(new Complex[][]{g, f}));
        assertArrayEquals(expected, sol.getArray());
    }
}
