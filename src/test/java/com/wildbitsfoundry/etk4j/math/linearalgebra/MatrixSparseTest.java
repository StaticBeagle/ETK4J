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
        MatrixSparse sparseCSC = MatrixSparse.from2dArray(matrix, ConstantsETK.DOUBLE_EPS);
        assertEquals(1, sparseCSC.get(0, 0), 1e-12);
    }

    @Test
    public void testUnsafeGet() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2dArray(matrix, ConstantsETK.DOUBLE_EPS);
        sparseCSC.set(0, 0, 100);
        assertEquals(100, sparseCSC.unsafeGet(0, 0), 1e-12);
    }

    @Test
    public void testLUDecomposition() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2dArray(matrix, ConstantsETK.DOUBLE_EPS);
        double[] b = {1, 2, 0};
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
        MatrixSparse sparseCSC = MatrixSparse.from2dArray(matrix, ConstantsETK.DOUBLE_EPS);
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
        MatrixSparse sparseCSC = MatrixSparse.from2dArray(matrix, ConstantsETK.DOUBLE_EPS);
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
        MatrixSparse sparseCSC = MatrixSparse.from2dArray(matrix, ConstantsETK.DOUBLE_EPS);
        boolean isFull = sparseCSC.isFull();
        assertTrue(isFull);
    }
}
