package com.wildbitsfoundry.etk4j.math.linearalgebra;

import java.util.Random;

public final class Matrices {
    private Matrices() {
    }

    public static Matrix Vandermonde(double[] x, int rows, int cols) {
        double[][] V = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                V[i][j] = Math.pow(x[i], j);
            }
        }
        return new Matrix(V);
    }

    public static Matrix Identity(int rows, int cols) {
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

    public static Matrix Identity(int n) {
        return Matrices.Identity(n, n);
    }

    public static Matrix Random(int rows, int cols) {
        Random rand = new Random();
        double[] data = new double[rows * cols];

        for(int i = 0; i < data.length; ++i) {
            data[i] = rand.nextDouble() * 100.0;
        }
        return new Matrix(data, rows, cols);
    }

    public static Matrix Companion(double[] coefs, int n) {
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

    public static Matrix Magic(int n) {
        Matrix magicMatrix = null;
        if (n == 1) {
            magicMatrix = new Matrix(n, n);
            magicMatrix.set(0, 0, 1.0);
        } else if (n % 2 == 1) {    // n is odd
            magicMatrix = oddMagicMatrix(n);
        } else if (n % 4 == 0) {     // if n is even and double order
            magicMatrix = evenDoubleOrderMagicMatrix(n);
        } else {
            magicMatrix = evenSigleOrderMagicMatrix(n);
        }
        return magicMatrix;
    }

    private static Matrix oddMagicMatrix(int n) {
        Matrix magicMatrix = new Matrix(n, n);
        int i = 0, j = n / 2, n2 = n * n;

        for (int k = 1; k <= n2; ++k) {
            magicMatrix.set(i, j, k);

            i--;
            j++;
            if (k % n == 0) {
                i += 2;
                --j;
            } else if (j == n) {
                j -= n;
            } else if (i < 0) {
                i += n;
            }
        }
        return magicMatrix;
    }

    private static Matrix evenDoubleOrderMagicMatrix(int n) {
        Matrix magicMatrix = new Matrix(n, n);
        int i, j;

        Matrix I = new Matrix(n, n);
        Matrix J = new Matrix(n, n);

        //prepare I, J
        int index = 1;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                I.set(i, j, (int) ((i + 1) % 4) / 2);
                J.set(j, i, (int) ((i + 1) % 4) / 2);
                magicMatrix.set(i, j, index);
                index++;
            }
        }

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (I.get(i, j) == J.get(i, j))
                    magicMatrix.set(i, j, n * n + 1 - magicMatrix.get(i, j));
            }
        }
        return magicMatrix;
    }

    private static Matrix evenSigleOrderMagicMatrix(int n) {

        int size = n * n;
        int halfN = n / 2;
        int subSquareSize = size / 4;

        Matrix subSquare = Magic(halfN);
        int[] quadrantFactors = {0, 2, 3, 1};
        double[][] result = new double[n][n];

        for (int r = 0; r < n; r++) {
            for (int c = 0; c < n; c++) {
                int quadrant = (r / halfN) * 2 + (c / halfN);
                result[r][c] = subSquare.get(r % halfN, c % halfN);
                result[r][c] += quadrantFactors[quadrant] * subSquareSize;
            }
        }

        int nColsLeft = halfN / 2;
        int nColsRight = nColsLeft - 1;

        for (int r = 0; r < halfN; r++)
            for (int c = 0; c < n; c++) {
                if (c < nColsLeft || c >= n - nColsRight
                        || (c == nColsLeft && r == nColsLeft)) {

                    if (c == 0 && r == nColsLeft)
                        continue;

                    double tmp = result[r][c];
                    result[r][c] = result[r + halfN][c];
                    result[r + halfN][c] = tmp;
                }
            }
        return new Matrix(result);
    }
}
