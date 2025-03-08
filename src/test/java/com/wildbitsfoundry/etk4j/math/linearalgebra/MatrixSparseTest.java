package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.constant.ConstantsETK;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;
import org.junit.Test;

import static org.junit.Assert.*;

public class MatrixSparseTest {

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
    public void testFrom2DArraySquareMatrix() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        double[] expected = {1, 4, 7, 2, 5, 8, 3, 6, 10};
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        assertArrayEquals(expected, sparseCSC.getArrayDense(), 1e-12);
    }

    @Test
    public void testFrom2DArrayTallMatrix() {
        double[][] matrix = {
                {1, 4},
                {2, 5},
                {3, 6},
        };
        double[] expected = {1, 4, 2, 5, 3, 6};
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        assertArrayEquals(expected, sparseCSC.getArrayDense(), 1e-12);
    }

    @Test
    public void testFrom2DArrayWideMatrix() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8}
        };
        double[] expected = {1, 4, 7, 2, 5, 8};
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        assertArrayEquals(expected, sparseCSC.getArrayDense(), 1e-12);
    }

    @Test
    public void testFrom2DArrayIdentityMatrix() {
        double[][] matrix = {
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1},
        };
        double[] expected = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        assertArrayEquals(expected, sparseCSC.getArrayDense(), 1e-12);
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
    public void testGetArrayDense() {
        MatrixSparse identity = MatrixSparse.Factory.identity(4);
        double[] expected = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
        assertArrayEquals(expected, identity.getArrayDense(), 1e-12);
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
    public void testSolveArrayRHS() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);
        double[] b = {1, 2, 0};
        MatrixSparse solution = sparseCSC.solve(b);
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
        MatrixSparse A = qrDecompositionSparse.getQ();

        double[] expected = {-0.2672612419124243, 0.8728715609439696, -0.4082482904638629, -0.5345224838248488,
                0.2182178902359924, 0.816496580927726, -0.8017837257372732, -0.4364357804719846, -0.4082482904638631};
        assertArrayEquals(expected, A.getArrayDense(), 1e-12);
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
        assertArrayEquals(expected, A.getArrayDense(), 1e-12);
    }

    @Test
    public void testGetH() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);

        // Perform Householder QR decomposition
        QRDecompositionSparse qrDecompositionSparse = new QRDecompositionSparse(sparseCSC);
        MatrixSparse A = qrDecompositionSparse.getH();

        double[] expected = {1.0, 0.0, 0.0, 0.42179344411906794, 1.0, 0.0, 0.632690166178602, 0.8597677535291113, 1.0};
        assertArrayEquals(expected, A.getArrayDense(), 1e-12);
    }

    @Test
    public void testQRDecompositionGetOriginalMatrixBack() {
        double[][] matrix = {
                {1, 4, 7},
                {2, 5, 8},
                {3, 6, 10},
        };
        MatrixSparse sparseCSC = MatrixSparse.from2DArray(matrix, ConstantsETK.DOUBLE_EPS);

        double[] expected = DoubleArrays.flatten(matrix);

        QRDecompositionSparse qrDecompositionSparse = new QRDecompositionSparse(sparseCSC);
        MatrixSparse A = qrDecompositionSparse.getQ().multiply(qrDecompositionSparse.getR());

        assertArrayEquals(expected, A.getArrayDense(), 1e-12);
    }

    @Test
    public void testQRDecompositionSolveSparse() {
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
    public void testQRDecompositionSolveSparseMatrixRHS() {
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

    @Test
    public void testGetL() {
        double[][] matrix = {
                {0.9767407521087846, 0.2794308404798361, 0.1685348957730604, 0.8646965674039528},
                {0.7198733467511716, 0.01506799866339292, 0.750545968546426, 0.8173900466200722},
                {0.3444764277139777, 0.8429046576872399, 0.6191817846917982, 0.1931980127831494},
                {0.2529620785541479, 0.7385598627121508, 0.3669833799842998, 0.4630012679866645}
        };

        MatrixSparse Arg = MatrixSparse.from2DArray(matrix);

        // Making it symmetric
        Arg = Arg.multiply(Arg.transpose());
        CholeskyDecompositionSparse choleskyDecompositionSparse = new CholeskyDecompositionSparse(Arg);

        double[] expected = {1.344696343497191, 1.145706730916891, 0.6272115798783218, 0.6809417540994949,
                0.661246176467614, 0.249108860226036, 0.1012637242504385, 0.8912991841476854, 0.6440431410415371, 0.2641588111155312};
        assertArrayEquals(expected, choleskyDecompositionSparse.getL().getNonZeroValues(), 1e-12);
    }

    @Test
    public void testCholeskyDecompositionGetOriginalMatrixBack() {
        double[][] matrix = {
                {0.9767407521087846, 0.2794308404798361, 0.1685348957730604, 0.8646965674039528},
                {0.7198733467511716, 0.01506799866339292, 0.750545968546426, 0.8173900466200722},
                {0.3444764277139777, 0.8429046576872399, 0.6191817846917982, 0.1931980127831494},
                {0.2529620785541479, 0.7385598627121508, 0.3669833799842998, 0.4630012679866645}
        };

        MatrixSparse Arg = MatrixSparse.from2DArray(matrix);

        // Making it symmetric
        Arg = Arg.multiply(Arg.transpose());

        CholeskyDecompositionSparse choleskyDecompositionSparse = new CholeskyDecompositionSparse(Arg);
        MatrixSparse L = choleskyDecompositionSparse.getL();
        MatrixSparse A = L.multiply(L.transpose());
        assertArrayEquals(Arg.getNonZeroValues(), A.getNonZeroValues(), 1e-12);
    }

    @Test
    public void tesCholeskyDecompositiontIsSPD() {
        double[][] matrix = {
                {0.9767407521087846, 0.2794308404798361, 0.1685348957730604, 0.8646965674039528},
                {0.7198733467511716, 0.01506799866339292, 0.750545968546426, 0.8173900466200722},
                {0.3444764277139777, 0.8429046576872399, 0.6191817846917982, 0.1931980127831494},
                {0.2529620785541479, 0.7385598627121508, 0.3669833799842998, 0.4630012679866645}
        };

        MatrixSparse Arg = MatrixSparse.from2DArray(matrix);

        // Making it symmetric
        Arg = Arg.multiply(Arg.transpose());
        CholeskyDecompositionSparse choleskyDecompositionSparse = new CholeskyDecompositionSparse(Arg);

        assertTrue(choleskyDecompositionSparse.isSPD());
    }

    @Test
    public void testCholeskyDecompositionSolve() {
        double[][] matrix = {
                {0.9767407521087846, 0.2794308404798361, 0.1685348957730604, 0.8646965674039528},
                {0.7198733467511716, 0.01506799866339292, 0.750545968546426, 0.8173900466200722},
                {0.3444764277139777, 0.8429046576872399, 0.6191817846917982, 0.1931980127831494},
                {0.2529620785541479, 0.7385598627121508, 0.3669833799842998, 0.4630012679866645}
        };

        MatrixSparse Arg = MatrixSparse.from2DArray(matrix);

        // Making it symmetric
        Arg = Arg.multiply(Arg.transpose());
        CholeskyDecompositionSparse choleskyDecompositionSparse = new CholeskyDecompositionSparse(Arg);

        MatrixSparse B = MatrixSparse.from2DArray(new double[][]{
                {90},
                {80},
                {60},
                {40}
        });

        MatrixSparse X = choleskyDecompositionSparse.solve(B);

        double[] expected = {110.8752235125764, -31.42490731167088, 209.5953908930501, -260.8458452894825};
        assertArrayEquals(expected, X.getNonZeroValues(), 1e-12);
    }

    @Test
    public void testCholeskyDecompositionSolveUseArray() {
        double[][] matrix = {
                {0.9767407521087846, 0.2794308404798361, 0.1685348957730604, 0.8646965674039528},
                {0.7198733467511716, 0.01506799866339292, 0.750545968546426, 0.8173900466200722},
                {0.3444764277139777, 0.8429046576872399, 0.6191817846917982, 0.1931980127831494},
                {0.2529620785541479, 0.7385598627121508, 0.3669833799842998, 0.4630012679866645}
        };

        MatrixSparse Arg = MatrixSparse.from2DArray(matrix);

        // Making it symmetric
        Arg = Arg.multiply(Arg.transpose());
        CholeskyDecompositionSparse choleskyDecompositionSparse = new CholeskyDecompositionSparse(Arg);

        double[] b = {90, 80, 60, 40};

        MatrixSparse X = choleskyDecompositionSparse.solve(b);

        double[] expected = {110.8752235125764, -31.42490731167088, 209.5953908930501, -260.8458452894825};
        assertArrayEquals(expected, X.getNonZeroValues(), 1e-12);
    }
}
