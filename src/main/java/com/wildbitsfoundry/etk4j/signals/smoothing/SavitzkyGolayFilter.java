package com.wildbitsfoundry.etk4j.signals.smoothing;

import com.wildbitsfoundry.etk4j.math.linearalgebra.LeastSquaresSolver;
import com.wildbitsfoundry.etk4j.math.linearalgebra.MatrixDense;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.Arrays;

// TODO write tests
public class SavitzkyGolayFilter {

    public enum Mode {
        INTERPOLATION,
        MIRROR,
        NEAREST,
        CONSTANT,
        WRAP
    }

    /**
     * Apply the Savitzky-Golay filter to smooth the input data.
     *
     * @param data       Input data array (noisy signal).
     * @param windowSize Size of the sliding window (must be odd).
     * @param polyOrder  Degree of the polynomial (must be less than window size).
     * @return Smoothed data array.
     */
    public static double[] applyFilter(double[] data, int windowSize, int polyOrder, Mode mode, double constant) {
        if (windowSize % 2 == 0) {
            throw new IllegalArgumentException("Window size must be odd.");
        }
        if (polyOrder >= windowSize) {
            throw new IllegalArgumentException("Polynomial order must be less than window size.");
        }

        int halfWindow = windowSize / 2;

        // Calculate coefficients
        double[] x = DoubleArrays.linSteps(-halfWindow, halfWindow, 1);

        MatrixDense A = MatrixDense.Factory.vandermonde(x, windowSize, polyOrder + 1);
        double[] coefficients = LeastSquaresSolver.solve(A, MatrixDense.Factory.identity(A.getRowCount())).getRow(0);

        switch (mode) {
            case INTERPOLATION:
                double[] smoothed =  convolve1d(data, coefficients, windowSize, Mode.MIRROR, constant);
                fitEdgesPolyFit(smoothed, windowSize, polyOrder, data);
                return smoothed;
            default:
                return convolve1d(data, coefficients, windowSize, mode, constant);
        }
    }

    private static void fitEdgesPolyFit(double[] smoothed, int windowSize, int polyOrder, double[] data) {
        int halfWindow = windowSize / 2;
        fitEdge(smoothed, 0, windowSize, 0 , halfWindow, polyOrder, data);

        int n = data.length;
        fitEdge(smoothed, n - windowSize, n, n - halfWindow, n, polyOrder, data);
    }

    private static void fitEdge(double[] smoothed, int windowStart, int windowStop, int x0, int xn, int polyOrder, double[] data) {
        double[] x = DoubleArrays.linSteps(0, windowStop - windowStart - 1);
        double[] y = new double[windowStop - windowStart];
        System.arraycopy(data, windowStart, y, 0, windowStop - windowStart);

        double[] polyCoefficients = Polynomial.polyFit(x, y, polyOrder).getCoefficients();
        double[] values = Polynomial.polyVal(polyCoefficients, DoubleArrays.linSteps(x0 - windowStart, xn - windowStart - 1));
        System.arraycopy(values, 0, smoothed, x0, values.length);
    }

    private static double[] convolve1d(double[] data, double[] coefficients, int windowSize, Mode mode, double constant) {
        switch (mode) {
            case MIRROR:
                return convolve1dMirror(data, coefficients, windowSize);
            case NEAREST:
                return convolve1dNearest(data, coefficients, windowSize);
            case CONSTANT:
                return convolve1dConstant(data, coefficients, windowSize, constant);
            case WRAP:
                return convolve1dWrap(data, coefficients, windowSize);
            default:
                throw new IllegalArgumentException("Mode");
        }
    }

    private static double[] convolve1dMirror(double[] data, double[] coefficients, int windowSize) {
        int halfWindow = windowSize / 2;
        double[] smoothed = new double[data.length];
        // Apply the filter
        for (int i = 0; i < data.length; i++) {
            double result = 0.0;
            for (int j = -halfWindow; j <= halfWindow; j++) {
                int index = i + j;

                // Handle boundary conditions (mirror edges)
                if (index < 0) {
                    index = -index; // Reflect at the start
                } else if (index >= data.length) {
                    index = 2 * data.length - index - 2; // Reflect at the end
                }
                result += coefficients[j + halfWindow] * data[index];
            }
            smoothed[i] = result;
        }
        return smoothed;
    }

    private static double[] convolve1dNearest(double[] data, double[] coefficients, int windowSize) {
        int halfWindow = windowSize / 2;
        double[] smoothed = new double[data.length];
        // Apply the filter
        for (int i = 0; i < data.length; i++) {
            double result = 0.0;
            for (int j = -halfWindow; j <= halfWindow; j++) {
                int index = i + j;

                // Handle boundary conditions
                if (index < 0) {
                    result += coefficients[j + halfWindow] * data[0];
                } else if (index >= data.length) {
                    result += coefficients[j + halfWindow] * data[data.length - 1];
                } else {
                    result += coefficients[j + halfWindow] * data[index];
                }
            }
            smoothed[i] = result;
        }
        return smoothed;
    }

    private static double[] convolve1dConstant(double[] data, double[] coefficients, int windowSize, double constant) {
        int halfWindow = windowSize / 2;
        double[] smoothed = new double[data.length];
        // Apply the filter
        for (int i = 0; i < data.length; i++) {
            double result = 0.0;
            for (int j = -halfWindow; j <= halfWindow; j++) {
                int index = i + j;

                // Handle boundary conditions
                if (index < 0) {
                    result += coefficients[j + halfWindow] * constant;
                } else if (index >= data.length) {
                    result += coefficients[j + halfWindow] * constant;
                } else {
                    result += coefficients[j + halfWindow] * data[index];
                }
            }
            smoothed[i] = result;
        }
        return smoothed;
    }

    private static double[] convolve1dWrap(double[] data, double[] coefficients, int windowSize) {
        int halfWindow = windowSize / 2;
        double[] smoothed = new double[data.length];
        // Apply the filter
        for (int i = 0; i < data.length; i++) {
            double result = 0.0;
            for (int j = -halfWindow; j <= halfWindow; j++) {
                int index = i + j;

                // Handle boundary conditions
                if (index < 0) {
                    index = data.length + index;
                } else if (index >= data.length) {
                    index = index - data.length;
                }
                result += coefficients[j + halfWindow] * data[index];
            }
            smoothed[i] = result;
        }
        return smoothed;
    }


    public static void main(String[] args) {
// Example noisy data
        double[] data = {2, 2, 5, 2, 1, 0, 1, 4, 9};
        int windowSize = 5; // Odd window size
        int polyOrder = 2; // Polynomial degree

        double[] smoothed = applyFilter(data, windowSize, polyOrder, Mode.INTERPOLATION, 5);

        System.out.println("Original data: " + Arrays.toString(data));
        System.out.println("Smoothed data: " + Arrays.toString(smoothed));
    }
}
