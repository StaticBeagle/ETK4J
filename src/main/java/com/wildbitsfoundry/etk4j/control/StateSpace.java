package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.linearalgebra.EigenvalueDecomposition;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.math.linearalgebra.NonSquareMatrixException;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.NumArrays;

import java.util.Arrays;

public class StateSpace extends LinearTimeInvariantSystem {

    private Matrix A;
    private Matrix B;
    private Matrix C;
    private Matrix D;

    public StateSpace() {
        this.A = Matrices.empty();
        this.B = Matrices.empty();
        this.C = Matrices.empty();
        this.D = Matrices.empty();
    }

    public StateSpace(double[][] A, double[][] B, double[][] C, double[][] D) {
        this.A = new Matrix(A);
        this.B = new Matrix(B);
        this.C = new Matrix(C);
        this.D = new Matrix(D);
    }

    public StateSpace(Matrix A, Matrix B, Matrix C, Matrix D) {
        this.A = new Matrix(A);
        this.B = new Matrix(B);
        this.C = new Matrix(C);
        this.D = new Matrix(D);
    }

    public Matrix getA() {
        return new Matrix(A);
    }

    public Matrix getB() {
        return new Matrix(B);
    }

    public Matrix getC() {
        return new Matrix(C);
    }

    public Matrix getD() {
        return new Matrix(D);
    }

    @Override
    public String toString() {
        return String.format("A:%n%s%nB:%n%s%nC:%n%s%nD%n%s%n", A, B, C, D);
    }

    @Override
    public StateSpace toStateSpace() {
        return this;
    }

    /**
     * https://github.com/scipy/scipy/blob/47bb6febaa10658c72962b9615d5d5aa2513fa3a/scipy/signal/lti_conversion.py#L196
     * @return
     */
    public TransferFunction[] toTransferFunction(int input) {
        //TODO normalizeABCD
        int nin = D.getColumnCount();
        int nout = D.getRowCount();
        if(input >= nin) {
            throw new IllegalArgumentException("System does not have the input specified.");
        }
        // make SIMO from possibly MIMO system.
        B = B.subMatrix(0, B.getRowCount() - 1, input, input);
        D = D.subMatrix(0, D.getRowCount() - 1, input, input);

        double[][] num;
        double[] den;

        if(!A.isSquared()) {
            throw new NonSquareMatrixException("Matrix A must be a square Matrix.");
        }
        den = A.poly();

        if(B.getRowCount() * B.getColumnCount() == 0 && C.getRowCount() * C.getColumnCount() == 0) {
            if(D.getRowCount() * D.getColumnCount() == 0 && A.getRowCount() * A.getColumnCount() == 0) {
                den = new double[0];
            }
            return new TransferFunction[] { new TransferFunction(D.getArray(), den)};
        }

        int numStates = A.getRowCount();
        num = new double[nout][numStates + 1];
        TransferFunction[] tfs = new TransferFunction[nout];
        for(int k = 0; k < nout; ++k) {
            double[] Ck = C.getRow(k);
            double Dk = D.get(k, 0);
            num[k] = A.subtract(new Matrix(dot(B.getArray(), Ck))).poly();
            NumArrays.addElementWiseInPlace(num[k], NumArrays.multiplyElementWise(den, Dk - 1));
            tfs[k] = new TransferFunction(num[k], den);
        }
        return tfs;
    }

    private static double[][] dot(double[] a, double[] b) {
        double[][] result = new double[a.length][b.length];
        for(int i = 0; i < a.length; ++i) {
            result[i] = NumArrays.multiplyElementWise(b, a[i]);
        }
        return result;
    }

    @Override
    public TransferFunction toTransferFunction() {
        return this.toTransferFunction(0)[0];
    }

    @Override
    public ZeroPoleGain toZeroPoleGain() {
        return this.toTransferFunction().toZeroPoleGain();
    }

    public TimeResponse simulateTimeResponse(double[][] input, double[] time) {
        return simulateTimeResponse(input, time, IntegrationMethod.INTERPOLATION);
    }

    public TimeResponse simulateTimeResponse(double[][] input, double[] time, double[] initialConditions) {
        return simulateTimeResponse(input, time, initialConditions, IntegrationMethod.INTERPOLATION);
    }

    public TimeResponse simulateTimeResponse(double[][] input, double[] time, IntegrationMethod integrationMethod) {
        return simulateTimeResponse(input, time, null, IntegrationMethod.INTERPOLATION);
    }

    public TimeResponse simulateTimeResponse(double[][] input, double[] time, double[] initialConditions,
                                             IntegrationMethod integrationMethod) {
        return lsim(input, time, initialConditions, this, integrationMethod);
    }

    public static void main(String[] args) {
        double[][] A = {{-2, -1}, {1, 0}};
        double[][] B = {{1}, {0}};
        double[][] C = {{1, 2}};
        double[][] D = {{1}};

        StateSpace ss = new StateSpace(A, B, C, D);
        TransferFunction tf = ss.toTransferFunction();
        System.out.println(tf);

        A = new double[][]{{-1}};
        B = new double[][]{{1}};
        C = new double[][]{{1}};
        D = new double[][]{{0}};

        ss = new StateSpace(A, B, C, D);
        tf = ss.toTransferFunction();
        System.out.println(tf);

        A = new double[][] {{-2, -1, 3}, {1, 0, 5}, {4, 5, 10}};
        B = new double[][]{{1, 2, 4}, {0, 6, 7}, {9, 10, 22}};
        C = new double[][]{{1, 2, 0}, {0, 1, 0}, {0, 0, 1}};
        D = new double[][]{{0, 0, 1}, {2, 3, 4}, {5, 6, 7}};

        ss = new StateSpace(A, B, C, D);
        TransferFunction[] tfs = ss.toTransferFunction(0);
        Arrays.stream(tfs).forEach(System.out::println);
    }
}
