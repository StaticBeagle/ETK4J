package com.wildbitsfoundry.etk4j.math.linearalgebra;

import java.util.Random;

// TODO move this to Matrix, is this class even necessary
public final class Matrices {
    private Matrices() {
    }

    public static Matrix vandermonde(double[] x, int rows, int cols) {
        double[][] V = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                V[i][j] = Math.pow(x[i], j);
            }
        }
        return new Matrix(V);
    }

    public static Matrix identity(int rows, int cols) {
        double[] data = new double[rows * cols];
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (i == j) {
                    data[i * cols + j] = 1.0;
                }
            }
        }
        return new Matrix(data, rows, cols);
    }

    public static Matrix identity(int n) {
        return Matrices.identity(n, n);
    }

    public static Matrix random(int n) {
        return random(n, n);
    }

    public static Matrix random(int rows, int cols) {
        Random rand = new Random();
        double[] data = new double[rows * cols];

        for (int i = 0; i < data.length; ++i) {
            data[i] = rand.nextDouble() * 100.0;
        }
        return new Matrix(data, rows, cols);
    }

    public static Matrix companion(double[] coefs, int n) {
        // Construct the companion matrix
        Matrix c = new Matrix(n, n);

        double a = 1.0 / coefs[0];
        for (int i = 0; i < n; i++) {
            c.set(0, n - 1 - i, -coefs[n - i] * a);
        }
        for (int i = 1; i < n; i++) {
            c.set(i, i - 1, 1);
        }
        return c;
    }

    public static Matrix empty() {
        return new Matrix(0, 0);
    }

    public static Matrix magic(int n) {
        Matrix magicMatrix = null;
        if (n == 1) {
            magicMatrix = new Matrix(n, n);
            magicMatrix.set(0, 0, 1.0);
        } else if (n == 2) {
            return Matrices.empty();
        }
        double[][] M = new double[n][n];

        // Odd order
        if ((n % 2) == 1) {
            int a = (n + 1) / 2;
            int b = (n + 1);
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < n; i++) {
                    M[i][j] = n * ((i + j + a) % n) + ((i + 2 * j + b) % n) + 1;
                }
            }
            // Doubly Even Order
        } else if ((n % 4) == 0) {
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < n; i++) {
                    if (((i + 1) / 2) % 2 == ((j + 1) / 2) % 2) {
                        M[i][j] = n * n - n * i - j;
                    } else {
                        M[i][j] = n * i + j + 1;
                    }
                }
            }
            // Singly Even Order
        } else {
            int p = n / 2;
            int k = (n - 2) / 4;
            Matrix A = magic(p);
            for (int j = 0; j < p; j++) {
                for (int i = 0; i < p; i++) {
                    double aij = A.get(i, j);
                    M[i][j] = aij;
                    M[i][j + p] = aij + 2 * p * p;
                    M[i + p][j] = aij + 3 * p * p;
                    M[i + p][j + p] = aij + p * p;
                }
            }
            for (int i = 0; i < p; i++) {
                for (int j = 0; j < k; j++) {
                    double t = M[i][j];
                    M[i][j] = M[i + p][j];
                    M[i + p][j] = t;
                }
                for (int j = n - k + 1; j < n; j++) {
                    double t = M[i][j];
                    M[i][j] = M[i + p][j];
                    M[i + p][j] = t;
                }
            }
            double t = M[k][0];
            M[k][0] = M[k + p][0];
            M[k + p][0] = t;

            t = M[k][k];
            M[k][k] = M[k + p][k];
            M[k + p][k] = t;
        }
        return new Matrix(M);
    }
}
