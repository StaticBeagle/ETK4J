package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import org.junit.Test;

import static org.junit.Assert.*;

public class MatrixSparseTest {

    // TODO write tests for constructors

    @Test
    public void testGet() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        assertEquals(1, sparseCSC.get(0, 0), 1e-12);
    }

    @Test
    public void testUnsafeGet() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        sparseCSC.set(0, 0, 100);
        assertEquals(100, sparseCSC.unsafeGet(0, 0), 1e-12);
    }

    @Test
    public void testAdd() {
        double[][] matrixA = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };

        double[][] matrixB = {
                {5, 9, 13},
                {3, 6, 7},
                {1, 8, 10},
        };

        MatrixSparse A = MatrixSparse.from2DArray(matrixA, ConstantsETK.DOUBLE_EPS);
        MatrixSparse B = MatrixSparse.from2DArray(matrixB, ConstantsETK.DOUBLE_EPS);
        MatrixDense C = A.add(B).toDense();

        double[] expected = {6, 13, 20, 5, 11, 15, 4, 14, 20};
        assertArrayEquals(expected, C.getArray(), 1e-12);
    }

    @Test
    public void testSubtract() {
        double[][] matrixA = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };

        double[][] matrixB = {
                {5, 9, 13},
                {3, 6, 7},
                {1, 8, 10},
        };

        MatrixSparse A = MatrixSparse.from2DArray(matrixA, ConstantsETK.DOUBLE_EPS);
        MatrixSparse B = MatrixSparse.from2DArray(matrixB, ConstantsETK.DOUBLE_EPS);
        MatrixDense C = A.subtract(B).toDense();

        double[] expected = {-4, -5, -6, -1, -1, 1, 2, -2, 0};
        assertArrayEquals(expected, C.getArray(), 1e-12);
    }

    @Test
    public void testMultiply() {
        double[][] matrixA = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };

        double[][] matrixB = {
                {5, 9, 13},
                {3, 6, 7},
                {1, 8, 10},
        };

        MatrixSparse A = MatrixSparse.from2DArray(matrixA, ConstantsETK.DOUBLE_EPS);
        MatrixSparse B = MatrixSparse.from2DArray(matrixB, ConstantsETK.DOUBLE_EPS);
        MatrixDense C = A.multiply(B).toDense();

        double[] expected = {24, 89, 111, 33, 112, 141, 43, 143, 181};
        assertArrayEquals(expected, C.getArray(), 1e-12);
    }

    @Test
    public void testMultiplyScalar() {
        double[][] matrixA = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };

        MatrixSparse A = MatrixSparse.from2DArray(matrixA, ConstantsETK.DOUBLE_EPS);
        MatrixDense C = A.multiply(2).toDense();

        double[] expected = {2, 8, 14, 4, 10, 16, 6, 12, 20};
        assertArrayEquals(expected, C.getArray(), 1e-12);
    }

    @Test
    public void testLUDecomposition() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        double[] b = {1, 2, 0};
        LUDecompositionSparse sparseLU = new LUDecompositionSparse(sparseCSC);
        MatrixSparse solution = sparseLU.solve(b);
        assertArrayEquals(new double[]{-2, 6.000000000000005, -3.0000000000000027}, new double[]{solution.get(0, 0),
                solution.get(1, 0), solution.get(2, 0)}, 1e-12);
    }

    @Test
    public void testLUDecompositionSparseMatrixRHS() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        MatrixSparse b = MatrixSparse.from2DArray(new double[][]{{1}, {2}, {ConstantsETK.DOUBLE_EPS}}, 0);
        LUDecompositionSparse sparseLU = new LUDecompositionSparse(sparseCSC);
        MatrixSparse solution = sparseLU.solve(b);
        assertArrayEquals(new double[]{-2, 6.000000000000005, -3.0000000000000027}, new double[]{solution.get(0, 0),
                solution.get(1, 0), solution.get(2, 0)}, 1e-12);
    }

    @Test
    public void unsafeSet() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        sparseCSC.unsafeSet(2, 2, 100);
        assertEquals(100, sparseCSC.unsafeGet(2, 2), 1e-12);
    }

    @Test
    public void testReshape() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        sparseCSC.reshape(10, 10);
        assertArrayEquals(new double[]{10, 10}, new double[]{sparseCSC.getRowCount(), sparseCSC.getColumnCount()}, 1e-12);
    }

    @Test
    public void testIsFull() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        boolean isFull = sparseCSC.isFull();
        assertTrue(isFull);
    }
}
