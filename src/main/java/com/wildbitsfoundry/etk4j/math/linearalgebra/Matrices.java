package com.wildbitsfoundry.etk4j.math.linearalgebra;

public final class Matrices {
    private Matrices() {
    }

    public static Matrix fwdSubsSolve(Matrix L, Matrix B) {
        final int nx = B.getColumnCount();
        final int m = L.getRowCount();
        final int n = L.getColumnCount();
        double[] t = L.getArray();
        double[] X = B.getArrayCopy();

        for (int j = 0; j < nx; ++j) {
            for (int i = 0; i < m; ++i) {
                for (int k = 0; k < i; ++k) {
                    X[i * nx + j] -= X[k * nx + j] * t[i * n + k];
                }
                X[i * nx + j] /= t[i * n + i];
            }
        }
        return new Matrix(X, m, nx);
    }

//	public static double[] backSubtitution(Matrix U, double[] b) {
//		final int length = b.length;
//		double[] x = new double[length];
//		for (int i = length - 1, index; i >= 0; --i) {
//			x[i] = b[i];
//			index = i * U._cols;
//			for (int j = i + 1; j < length; ++j) {
//				x[i] -= U._data[index + j] * x[j];
//			}
//			x[i] /= U._data[index + i];
//		}
//		return x;
//	}

    public static void tridiagonalLDLTSolve(double[] D, double[] L, double[] b, int length, int offset) {
        int n = length;
        for (int i = offset; i < n - 1; ++i) {
            double ui = L[i];
            L[i] /= D[i];
            D[i + 1] -= ui * L[i];
            b[i + 1] -= L[i] * b[i];
        }

        b[n - 1] /= D[n - 1];
        for (int i = n - 2; i >= offset; --i) {
            b[i] = b[i] / D[i] - L[i] * b[i + 1];
        }
    }

    public static void tridiagonalLDLTSolve(double[] D, double[] L, double[] b) {
        tridiagonalLDLTSolve(D, L, b, b.length, 0);
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

    private static Matrix oddMagicMatrix(int n) {
        Matrix magicMatrix = new Matrix(n, n);
        int i = 1, j = n / 2 + 1, n2 = n * n;

        for (int k = 1; k <= n2; ++k) {
            magicMatrix.set(i, j, k);

            i--;
            j++;
            if (k % n == 0) {
                i += 2;
                --j;
            } else {
                if (j == n + 1)
                    j -= n;
                else if (i < 1)
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
        for (i = 1; i <= n; i++) {
            for (j = 1; j <= n; j++) {
                I.set(i, j, (int) (i % 4) / 2);
                J.set(j, i, (int) (i % 4) / 2);
                magicMatrix.set(i, j, index);
                index++;
            }
        }

        for (i = 1; i <= n; i++) {
            for (j = 1; j <= n; j++) {
                if (I.get(i, j) == J.get(i, j))
                    magicMatrix.set(i, j, n * n + 1 - magicMatrix.get(i, j));
            }
        }
        return magicMatrix;
    }

    private static Matrix evenMagicMatrix(int n) {
        Matrix magicMatrix = new Matrix(n, n);
        int i, j, k;
        int w = n / 2;

        Matrix Mat = Magic(w);

        for (i = 1; i <= w; i++) {
            for (j = 1; j <= w; j++) {
                magicMatrix.set(i, j, Mat.get(i, j));
                magicMatrix.set(i + w, j, Mat.get(i, j) + 3 * w * w);
                magicMatrix.set(i, j + w, Mat.get(i, j) + 2 * w * w);
                magicMatrix.set(i + w, j + w, Mat.get(i, j) + w * w);
            }
        }
        if (n == 2) {
            return magicMatrix;
        }

        Matrix J = new Matrix(w, w - 1);
        k = (n - 2) / 4;
        for (i = 1; i <= w; i++) {
            for (j = 1; j <= k; j++) {
                J.set(i, j, magicMatrix.get(i, j));
                magicMatrix.set(i, j, magicMatrix.get(i + w, j));
                magicMatrix.set(i + w, j, J.get(i, j));

            }
        }

        for (i = 1; i <= w; i++) {
            for (j = 1; j < k; j++) {
                J.set(i, j, magicMatrix.get(i, n - j + 1));
                magicMatrix.set(i, n - j + 1, magicMatrix.get(i + w, n - j + 1));
                magicMatrix.set(i + w, n - j + 1, J.get(i, j));

            }
        }

        double temp;
        k = w / 2 + 1;
        temp = magicMatrix.get(k, 1);
        magicMatrix.set(k, 1, magicMatrix.get(k + w, 1));
        magicMatrix.set(k + w, 1, temp);

        temp = magicMatrix.get(k, k);
        magicMatrix.set(k, k, magicMatrix.get(k + w, k));
        magicMatrix.set(k + w, k, temp);

        return magicMatrix;
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
            magicMatrix = evenMagicMatrix(n);
        }
        return magicMatrix;
    }
}
