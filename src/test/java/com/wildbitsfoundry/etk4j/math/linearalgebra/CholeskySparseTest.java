package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;
// TODO add calculate R and L for all Cholesky decompositions
public class CholeskySparseTest {

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
    public void testGetOriginalMatrixBack() {
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
    public void testIsSPD() {
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
    public void testSolve() {
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

        // TODO add solve method for double[]
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
}
