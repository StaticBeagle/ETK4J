package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;

public class ComplexMatrixTest {

    @Test
    public void testInverse() {
        MatrixDense A = new MatrixDense(new double[][]{{-2, -1}, {1, 0}});

        double w = 100;
        Complex jw = Complex.fromImaginary(w);
        ComplexMatrixDense cm = ComplexMatrixDense.Factory.identity(A.getRowCount()).multiply(jw).subtract(A);
        ComplexMatrixDense inv = cm.inv();
        double[] real = {1.999600059992001E-4, 9.99700049993001E-5, -9.99700049993001E-5, 1.9996000599920008E-8};
        double[] imag = {-0.00999700049993001, 1.9996000599920007E-6, -1.9996000599920007E-6, -0.010000999700049992};

        Complex[] expected = ComplexArrays.zip(real, imag);
        assertArrayEquals(expected, inv.getArray());
    }

    @Test
    public void testQR() { // TODO move qr tests here?
        //ComplexMatrixDense A = new ComplexMatrixDense();
    }

    @Test
    public void testEigenValueDecomposition() {
        Complex[][] matrix = {
                {new Complex(65, 24), new Complex(35, 55), new Complex(40, 89), new Complex(69, 64)},
                {new Complex(99, 66), new Complex(64, 87), new Complex(37, 27), new Complex(2, 32)},
                {new Complex(39, 50), new Complex(48, 45), new Complex(35, 69), new Complex(90, 3)},
                {new Complex(30, 82), new Complex(93, 40), new Complex(87, 99), new Complex(17, 44)}
        };

        ComplexMatrixDense A = ComplexMatrixDense.from2DArray(matrix);
        ComplexEigenvalueDecompositionDense eig = new ComplexEigenvalueDecompositionDense(A);

        // Test eigenvalues
        Complex[] expected = {new Complex(215.25480299438541, 215.46132164576053), new Complex(), new Complex(), new Complex(),
                new Complex(), new Complex(17.957628223412968, -61.174250332479), new Complex(), new Complex(),
                new Complex(), new Complex(), new Complex(15.288613829263491, 68.94768884988993), new Complex(),
                new Complex(), new Complex(), new Complex(), new Complex(-67.50104504706171, 0.7652398368285932)

        };
        assertArrayEquals(expected, eig.getD().getArray());
    }

    @Test
    public void testVerifyEigenvectors() {
        Complex[][] matrix = {
                {new Complex(65, 24), new Complex(35, 55), new Complex(40, 89), new Complex(69, 64)},
                {new Complex(99, 66), new Complex(64, 87), new Complex(37, 27), new Complex(2, 32)},
                {new Complex(39, 50), new Complex(48, 45), new Complex(35, 69), new Complex(90, 3)},
                {new Complex(30, 82), new Complex(93, 40), new Complex(87, 99), new Complex(17, 44)}
        };

        ComplexMatrixDense A = ComplexMatrixDense.from2DArray(matrix);
        ComplexEigenvalueDecompositionDense eig = new ComplexEigenvalueDecompositionDense(A);

        // Verify that the results satisfy A*V = V*D
        ComplexMatrixDense result = A.multiply(eig.getV()).subtract(eig.getV().multiply(eig.getD()));

        String[] real = Arrays.stream(result.getArray()).map(re -> String.valueOf(re.real())).toArray(String[]::new);
        String[] imag = Arrays.stream(result.getArray()).map(re -> String.valueOf(re.imag())).toArray(String[]::new);

        String[] expectedReal = {"-3.552713678800501E-15", "6.039613253960852E-14", "-1.9539925233402755E-14",
                "-1.0436096431476471E-14", "-8.171241461241152E-14", "5.6843418860808015E-14", "0.0",
                "6.039613253960852E-14", "2.8421709430404007E-14", "4.440892098500626E-14", "2.8421709430404007E-14",
                "7.105427357601002E-15", "-1.0658141036401503E-14", "7.72715225139109E-14", "-7.105427357601002E-15", "0.0"};
        String[] expectedImag = {"-2.8421709430404007E-14", "2.1316282072803006E-14", "-8.881784197001252E-15",
                "-7.105427357601002E-14", "-1.7053025658242404E-13", "-8.881784197001252E-15", "-1.4210854715202004E-14",
                "-5.240252676230739E-14", "0.0", "-7.105427357601002E-15", "-3.197442310920451E-14", "-2.6645352591003757E-14",
                "-1.1368683772161603E-13", "-7.105427357601002E-15", "-2.8421709430404007E-14", "-3.197442310920451E-14"};

        // Test eigenvalues
        assertArrayEquals(expectedReal, real);
        assertArrayEquals(expectedImag, imag);
    }

    @Test
    public void testSingularValueDecompositionSingularValues() {
        Complex[][] matrix = {
                {new Complex(65, 24), new Complex(35, 55), new Complex(40, 89), new Complex(69, 64)},
                {new Complex(99, 66), new Complex(64, 87), new Complex(37, 27), new Complex(2, 32)},
                {new Complex(39, 50), new Complex(48, 45), new Complex(35, 69), new Complex(90, 3)},
                {new Complex(30, 82), new Complex(93, 40), new Complex(87, 99), new Complex(17, 44)}
        };
        ComplexMatrixDense A = new ComplexMatrixDense(matrix);
        ComplexSingularValueDecompositionDense svd = new ComplexSingularValueDecompositionDense(A);
        Complex[] expected = {Complex.fromReal(312.12634848312956), Complex.fromReal(106.49976593122118),
                Complex.fromReal(79.32524103615428), Complex.fromReal(35.10624693326205)};
        Complex[] actual = svd.S.diag();
        assertArrayEquals(expected, actual);
    }

    @Test
    public void testVerifySingularValueDecomposition() {
        Complex[][] matrix = {
                {new Complex(65, 24), new Complex(35, 55), new Complex(40, 89), new Complex(69, 64)},
                {new Complex(99, 66), new Complex(64, 87), new Complex(37, 27), new Complex(2, 32)},
                {new Complex(39, 50), new Complex(48, 45), new Complex(35, 69), new Complex(90, 3)},
                {new Complex(30, 82), new Complex(93, 40), new Complex(87, 99), new Complex(17, 44)}
        };
        ComplexMatrixDense A = new ComplexMatrixDense(matrix);
        ComplexSingularValueDecompositionDense svd = new ComplexSingularValueDecompositionDense(A);
        // U * S * VH
        Complex[] actual = svd.getU().multiply(svd.getS()).multiply(svd.getV().conjugateTranspose()).getArray();
        Complex[] expected = A.getArray();

        boolean isClose = true;
        for(int i = 0; i < expected.length; i++) {
            isClose &= MathETK.isClose(expected[i].real(), actual[i].real(), 1e-12, 0) &&
                    MathETK.isClose(expected[i].imag(), actual[i].imag(), 1e-12, 0);
        }
        assertTrue(isClose);
    }
}
