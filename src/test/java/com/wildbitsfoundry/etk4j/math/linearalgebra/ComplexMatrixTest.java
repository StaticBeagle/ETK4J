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
        double[] expected = {312.12634848312956, 106.49976593122118, 79.32524103615428, 35.10624693326205};
        double[] actual = svd.getS().diag();
        assertArrayEquals(expected, actual, 1e-12);
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
        ComplexMatrixDense S = ComplexMatrixDense.Factory.zeros(A.getRowCount(), A.getColumnCount());
        for (int i = 0; i < A.getRowCount(); i++) {
            S.unsafeSet(i, i, Complex.fromReal(svd.getS().diag()[i]));
        }
        // U * S * VH
        Complex[] actual = svd.getU().multiply(S).multiply(svd.getV().conjugateTranspose()).getArray();
        Complex[] expected = A.getArray();

        boolean isClose = true;
        for (int i = 0; i < expected.length; i++) {
            isClose &= MathETK.isClose(expected[i].real(), actual[i].real(), 1e-12, 0) &&
                    MathETK.isClose(expected[i].imag(), actual[i].imag(), 1e-12, 0);
        }
        assertTrue(isClose);
    }

    @Test
    public void testUnderDeterminedSingularValueDecomposition() {
        double[][] matrix = {
                {65, 35, 40, 69},
                {99, 64, 37, 2}
        };
        ComplexMatrixDense A = ComplexMatrixDense.fromRealMatrix(MatrixDense.from2DArray(matrix));
        ComplexSingularValueDecompositionDense svd = new ComplexSingularValueDecompositionDense(A);
        ComplexMatrixDense S = ComplexMatrixDense.Factory.zeros(A.getRowCount(), A.getColumnCount());
        for (int i = 0; i < A.getRowCount(); i++) {
            S.unsafeSet(i, i, Complex.fromReal(svd.getS().diag()[i]));
        }
        // U * S * VH
        Complex[] actual = svd.getU().multiply(S).multiply(svd.getV().conjugateTranspose()).getArray();
        Complex[] expected = A.getArray();

        boolean isClose = true;
        for (int i = 0; i < expected.length; i++) {
            isClose &= MathETK.isClose(expected[i].real(), actual[i].real(), 1e-12, 0) &&
                    MathETK.isClose(expected[i].imag(), actual[i].imag(), 1e-12, 0);
        }
        assertTrue(isClose);
    }

    @Test
    public void testPinv() {
        Complex[][] matrix = {
                {new Complex(65, 24), new Complex(35, 55), new Complex(40, 89), new Complex(69, 64)},
                {new Complex(99, 66), new Complex(64, 87), new Complex(37, 27), new Complex(2, 32)},
                {new Complex(39, 50), new Complex(48, 45), new Complex(35, 69), new Complex(90, 3)},
                {new Complex(30, 82), new Complex(93, 40), new Complex(87, 99), new Complex(17, 44)}
        };
        ComplexMatrixDense A = new ComplexMatrixDense(matrix).pinv();
        Complex[] expected = {new Complex(0.006606693008411897, 0.009028508690222058), new Complex(0.002241973628372189, -6.870241561510043E-4), new Complex(6.829853667463159E-4, -0.009483921279180879), new Complex(-0.005191153135118873, -0.0021323448816482748),
                new Complex(-0.011923334824432807, -0.0013132817812679708), new Complex(0.0030095244065256287, -0.008449216190560736), new Complex(0.00566933955476279, 0.005288234410094314), new Complex(0.005318529482423536, 0.0027846535437693197),
                new Complex(0.007083552797055823, -0.009010143220009651), new Complex(-3.468097572918551E-4, 0.007069615992465642), new Complex(-0.011046620188996124, 6.392546417408912E-4), new Complex(0.003990987778538757, -0.002040313883993867),
                new Complex(-0.0018998860419538865, -0.002784104098995203), new Complex(-0.0016404826783108246, -3.750688769298307E-4), new Complex(0.010151354796239605, 0.005957299401032975), new Complex(-0.0035784315548526803, -0.0024834440642732356)
        };
        assertArrayEquals(expected, A.getArray());
    }
}
