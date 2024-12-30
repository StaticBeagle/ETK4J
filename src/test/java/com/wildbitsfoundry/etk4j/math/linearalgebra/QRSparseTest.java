package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;
import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;

public class QRSparseTest {

    @Test
    public void testGetQ() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);

        // Perform Householder QR decomposition
        QRDecompositionSparse qrDecompositionSparse = new QRDecompositionSparse(sparseCSC);
        MatrixSparse A = qrDecompositionSparse.getQ(false);

        double[] expected = {-0.2672612419124243, 0.8728715609439696, -0.4082482904638629, -0.5345224838248488,
                0.2182178902359924, 0.816496580927726, -0.8017837257372732, -0.4364357804719846, -0.4082482904638631};
        double[] actual = new double[expected.length];

        for(int i = 0; i < A.getRowCount(); i++) {
            for(int j = 0; j < A.getColumnCount(); j++) {
                actual[A.getRowCount() * i + j] = A.unsafeGet(i, j);
            }
        }
        assertArrayEquals(expected, actual, 1e-12);
    }

    @Test
    public void testGetR() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);

        // Perform Householder QR decomposition
        QRDecompositionSparse qrDecompositionSparse = new QRDecompositionSparse(sparseCSC);
        MatrixSparse A = qrDecompositionSparse.getR();

        double[] expected = {-3.741657386773941, -8.552359741197581, -14.16484582135849, 0, 1.963961012123933,
                3.491486243775882, 0, 0, -0.4082482904638638};
        double[] actual = new double[expected.length];

        for(int i = 0; i < A.getRowCount(); i++) {
            for(int j = 0; j < A.getColumnCount(); j++) {
                actual[A.getRowCount() * i + j] = A.unsafeGet(i, j);
            }
        }
        assertArrayEquals(expected, actual, 1e-12);
    }

    @Test
    public void testGetV() {
        // TODO and decide on name V vs H
    }

    @Test
    public void testGetOriginalMatrixBack() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);

        double[] expected = DoubleArrays.flatten(matrix);
        double[] actual = new double[expected.length];

        QRDecompositionSparse qrDecompositionSparse = new QRDecompositionSparse(sparseCSC);
        MatrixSparse A = qrDecompositionSparse.getQ(false).multiply(qrDecompositionSparse.getR());
        for(int i = 0; i < A.getRowCount(); i++) {
            for(int j = 0; j < A.getColumnCount(); j++) {
                actual[A.getRowCount() * i + j] = A.unsafeGet(i, j);
            }
        }
        // Perform Householder QR decomposition
        assertArrayEquals(expected, actual, 1e-12);
    }

    @Test
    public void testSolveSparse() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        double[] b = {1, 2, 0};
        QRDecompositionSparse sparseQR = new QRDecompositionSparse(sparseCSC);
        MatrixSparse solution = sparseQR.solve(b);
        assertArrayEquals(new double[]{-2, 6.000000000000005, -3.0000000000000027}, new double[]{solution.get(0, 0),
                solution.get(1, 0), solution.get(2, 0)}, 1e-12);
    }

    @Test
    public void testSolveSparseMatrixRHS() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        MatrixSparse b = MatrixSparse.from2DArray(new double[][]{{1}, {2}, {ConstantsETK.DOUBLE_EPS}}, 0);
        QRDecompositionSparse sparseQR = new QRDecompositionSparse(sparseCSC);
        MatrixSparse solution = sparseQR.solve(b);
        assertArrayEquals(new double[]{-2, 6.000000000000005, -3.0000000000000027}, new double[]{solution.get(0, 0),
                solution.get(1, 0), solution.get(2, 0)}, 1e-12);
    }
}
