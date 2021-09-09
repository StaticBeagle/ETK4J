package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.math.linearalgebra.EigenvalueDecomposition;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.NumArrays;

import java.util.Arrays;

public abstract class LinearTimeInvariantSystem {

    protected abstract StateSpace toStateSpace();
    protected abstract TransferFunction toTransferFunction();
    protected abstract ZeroPoleGain toZeroPoleGain();

    protected TimeResponse lsim(double[][] input, double[] time, double[] initialConditions,
                              StateSpace ss) {
        double[][] U = NumArrays.transpose(input);

        if(U.length != time.length) {
            throw new IllegalArgumentException("The input array and the time array must have the same length.");
        }
        if(time.length == 0) {
            throw new IllegalArgumentException("The time array must have at least one element.");
        }

        Matrix A = ss.getA();
        Matrix B = ss.getB();
        Matrix C = ss.getC();
        Matrix D = ss.getD();

        final int noStates = A.getRowCount();
        final int noInputs = B.getColumnCount();
        final int noSteps = time.length;

        // initial conditions
        double[] x0 = initialConditions == null ? new double[noStates] : initialConditions;
        double[][] xOut = new double[noSteps][noStates];

        if(time[0] == 0.0) {
            xOut[0] = x0;
        } else if(time[0] > 0.0) {
            xOut[0] = NumArrays.dot(x0, A.transpose().multiply(time[0]).expm().getAs2DArray());
        } else {
            throw new IllegalArgumentException("Initial time must be non negative.");
        }

        if(noSteps > 1) {
            double dt = time[1] - time[0];
            double[] delta = new double[time.length - 2];
            for (int i = 1; i < time.length - 1; ++i) {
                delta[i - 1] = (time[i + 1] - time[i]) / dt;
            }
            if (!NumArrays.allClose(delta, 1.0)) {
                throw new NonUniformTimeStepsException();
            }

            A.multiplyEquals(dt);
            B.multiplyEquals(dt);
            double[][] M = new double[noStates + noInputs][];
            for (int i = 0; i < noStates; ++i) {
                M[i] = NumArrays.concatenate(A.getRow(i), B.getRow(i));
            }
            for(int i = noStates; i < noStates + noInputs; ++i) {
                M[i] = new double[noStates + noInputs];
            }

            Matrix expMT = new Matrix(M).transpose().expm();
            double[][] Ad = expMT.subMatrix(0, noStates - 1, 0, noStates - 1).getAs2DArray();
            double[][] Bd = expMT.subMatrix(noStates, expMT.getRowCount() - 1, 0, noStates - 1).getAs2DArray();
            for (int i = 1; i < noSteps; ++i) {
                xOut[i] = NumArrays.add(NumArrays.dot(xOut[i - 1], Ad), NumArrays.dot(U[i - 1], Bd));
            }
        }
        double[][] yOut = new double[noSteps][noStates];
        double[][] c = C.transpose().getAs2DArray();
        double[][] d = D.transpose().getAs2DArray();
        for(int i = 0; i < noSteps; ++i) {
            yOut[i] = NumArrays.dot(xOut[i], c);
            NumArrays.addElementWiseInPlace(yOut[i], NumArrays.dot(U[i], d));
        }
        return new TimeResponse(time, NumArrays.transpose(yOut), xOut);
    }

    public StepResponse step() {
        return step(100);
    }

    public StepResponse step(int numberOfPoints) {
        return stepResponse(null, null, numberOfPoints);
    }

    public StepResponse step(double[] initialConditions, int numberOfPoints) {
        return stepResponse(null, initialConditions, numberOfPoints);
    }

    public StepResponse step(double[] time) {
        return stepResponse(time, null, time.length);
    }

    public StepResponse step(double[] time, double[] initialConditions) {
        return stepResponse(time, initialConditions, time.length);
    }

    protected StepResponse stepResponse(double[] time, double[] initialConditions, int numberOfPoints) {
        StateSpace ss = this.toStateSpace();
        time = time == null ? defaultResponseTimes(ss.getA(), numberOfPoints) : time;
        double[][] U = new double[1][];
        U[0] = NumArrays.ones(time.length);
        TimeResponse lSim = lsim(U, time, initialConditions, ss);
        return new StepResponse(lSim.getTime(), lSim.getResponse()[0]);
    }

    protected double[] defaultResponseTimes(Matrix A, int numberOfPoints) {
        EigenvalueDecomposition eig = A.eig();
        double[] realEig = eig.getRealEigenvalues();
        for(int i = 0; i < realEig.length; ++i) {
            realEig[i] = Math.abs(realEig[i]);
        }
        double r = NumArrays.min(realEig);
        if(r == 0.0) {
            r = 1.0;
        }
        double tc = 1.0 / r;
        return NumArrays.linSpace(0.0, 7 * tc, numberOfPoints);
    }
}
