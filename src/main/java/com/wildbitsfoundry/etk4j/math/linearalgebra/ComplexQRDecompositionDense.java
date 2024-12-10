package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.Tuples;

import java.util.Arrays;

/***
 * References https://blogs.mathworks.com/cleve/2016/10/03/householder-reflections-and-the-qr-decomposition/
 */
public class ComplexQRDecompositionDense extends ComplexQRDecomposition<ComplexMatrixDense> {
    protected Complex[] _data;
    protected Complex[] _rdiag;

    public static void main(String[] args) {
// Define a 2x2 complex matrix
//        Complex[][] A = {
//                {new Complex(1, 1), new Complex(2, 2)},
//                {new Complex(3, 3), new Complex(4, 4)}
//        };
//
//// Perform Householder QR decomposition
//        ComplexQRDecompositionDense result = new ComplexQRDecompositionDense(new ComplexMatrixDense(A));
//
//// Print results
//        System.out.println("Matrix Q:");
//        System.out.println(result.getQ());
//
//        System.out.println("Matrix Q thin:");
//        System.out.println(result.getQThin());
//
//        System.out.println("Matrix R:");
//        System.out.println(result.getR());
//
//        System.out.println("Matrix H:");
//        System.out.println(result.getH());
//
//        System.out.println("Matrix Q * R:");
//        System.out.println(result.getQ().multiply(result.getR()));

        Complex[] x = {new Complex(1, 1), new Complex(2, 2)};
        Complex[] u = {new Complex(3, 3), new Complex(4, 4)};

        Complex[][] y = {x, u};
        //houseQR(new ComplexMatrixDense(y));

        Complex[] g = {new Complex(60, 70)};
        Complex[] f = {new Complex(80, 90)};

        ComplexMatrixDense sol = solve(new ComplexMatrixDense(y), new ComplexMatrixDense(new Complex[][]{g, f}));
        System.out.println(sol);
    }

    private static Complex sig(Complex u) {
        return u.sign().add(u.equals(new Complex()) ? 1 : 0);
    }

    private static Complex[] houseGen(Complex[] x) {
        double nu = ComplexArrays.norm(x);
        if(nu != 0) {
            Complex[] u = ComplexArrays.divideElementWise(x, nu);
            u[0] = u[0].add(sig(u[0]));
            return ComplexArrays.divideElementWise(u, Math.sqrt(u[0].abs()));
        } else {
            Complex[] u = ComplexArrays.deepCopy(x);
            u[1] = Complex.fromReal(Math.sqrt(2));
            return u;
        }
    }

    private static Complex[][] zeros(int m, int n) {
        Complex[][] zeros = new Complex[m][n];
        for(int i = 0; i < m; i++) {
            for(int j = 0; j < n; j++) {
                zeros[i][j] = new Complex();
            }
        }
        return zeros;
    }

    private static Complex[] zeros(int dim) {
        Complex[] zeros = new Complex[dim];
        for(int i = 0; i < dim; i++) {
            zeros[i] = new Complex();
        }
        return zeros;
    }

    private static Complex[][] calculateReflector(Complex[] u, Complex[][] x) {
        Complex[] uH = Arrays.stream(u).map(c -> c.conj()).toArray(Complex[]::new);
        Complex[] uHx = zeros(Math.max(x.length, x[0].length));
        for(int i = 0; i < x[0].length; i++) {
            for(int j = 0; j < x.length; j++) {
                // (uH * x)
                uHx[i].addEquals(uH[j].multiply(x[j][i]));
            }
        }
        Complex[][] uuHx = zeros(x.length, x[0].length);
        for(int i = 0; i < x.length; i++) {
            for(int j = 0; j < x[0].length; j++) {
                uuHx[i][j] = u[i].multiply(uHx[j]);
            }
        }

        Complex[][] H = new Complex[x.length][x[0].length];
        for(int i = 0; i < x.length; i++) {
            H[i] = ComplexArrays.subtractElementWise(x[i], uuHx[i]);
        }
        return H;
    }

    private static Complex[] getColumn(ComplexMatrixDense A, int col) {
        Complex[] result = new Complex[A.getRowCount()];
        for(int i = 0; i < result.length; i++) {
            result[i] = A.unsafeGet(i, col);
        }
        return result;
    }

    private static void setColumn(ComplexMatrixDense A, Complex[] values, int col, int startingRow) {
        for(int i = 0; i < A.getRowCount() - startingRow; i++) { // TODO could be endRow - startingRow if we pass a param
            A.set(i + startingRow, col, values[i]);
        }
    }

    private static void setColumn(ComplexMatrixDense A, Complex[] values, int col) {
        setColumn(A, values, col, 0);
    }

    private static void setColumn(ComplexMatrixDense A, Complex[] values, int row0, int row1, int col) {
        for(int i = row0; i < row1; i++) {
            A.set(i, col, values[i]);
        }
    }

//    private static ComplexMatrixDense getSubMatrix(ComplexMatrixDense A, int row0, int row1, int col0, int col1) {
//        ComplexMatrixDense result = new ComplexMatrixDense(row1 - row0, col1 - col0);
//        for(int i = row0; i < row1; i++) {
//            for(int j = col0; j < col1; j++) {
//                A.set
//            }
//        }
//    }

    private static Complex[][] getSubMatrix(ComplexMatrixDense A, int row0, int row1, int col0, int col1) {
        Complex[][] result = new Complex[row1 - row0][col1 - col0];
        for(int i = row0; i < row1; i++) {
            for(int j = col0; j < col1; j++) {
                result[i - row0][j - col0] = A.get(i, j);
            }
        }
        return result;
    }

    private static void setSubMatrix(ComplexMatrixDense A, Complex[][] values, int row0, int row1, int col0, int col1) {
//        if(col0 == col1) {
//            setColumn(A, values[0], row0, row1, col0);
//        } else {
//            for(int i = row0; i < row1; i++) {
//                for(int j = col0; j < col1; j++) {
//                    A.set(i, j, values[i - row0][j - col0]);
//                }
//            }
//        }
        for(int i = row0; i < row1; i++) {
            for(int j = col0; j < col1; j++) {
                A.set(i, j, values[i - row0][j - col0]);
            }
        }
        // TODO check for row0 == row1
    }

    public static Tuples.Tuple2<ComplexMatrixDense, ComplexMatrixDense> houseQR(ComplexMatrixDense A) {
        int m = A.getColumnCount();
        int n = A.getRowCount();
        ComplexMatrixDense R = A.copy();
        ComplexMatrixDense U = new ComplexMatrixDense(zeros(m, n));
        for(int i = 0; i < Math.min(m, n); i++) {
            Complex[] rCol = Arrays.copyOfRange(getColumn(R, i), i, m); // TODO get subColumn
            Complex[] u = houseGen(rCol);
            setColumn(U, u, i, i); // TODO set subColumn
            Complex[][] subR = getSubMatrix(R, i, m, i, n);
            Complex[][] H = calculateReflector(u, subR);
            setSubMatrix(R, H, i, m, i, n);
            for(int j = i + 1; j < m; j++) {
                R.unsafeSet(j, i, new Complex());
            }
        }
        System.out.println(R);
        System.out.println(U);

//        I = eye(size(U));
//        Q = house_apply(U,I)
        // TODO getQ
        ComplexMatrixDense I = ComplexMatrixDense.Factory.identity(U.getRowCount(), U.getColumnCount());
        ComplexMatrixDense Q = houseApply(U, I);

        // Check that Q is orthogonal
        System.out.println(Q.conjugateTranspose().multiply(Q));

        // Check Transpose
        System.out.println(houseApplyTranspose(U, I));

        // Check that we got A back
        System.out.println(Q.multiply(R));

        return new Tuples.Tuple2<>(U, R);
    }

    /*
     % Apply Householder reflections.
    % Z = house_apply(U,X), with U from house_qr
    % computes Q*X without actually computing Q.
     */
    public static ComplexMatrixDense houseApply(ComplexMatrixDense U, ComplexMatrixDense X) {
        ComplexMatrixDense Z = X.copy();
        int n = U.getColumnCount();
        Complex[][] H = new Complex[X.getRowCount()][X.getColumnCount()];
        // TODO temporary 2D array method
        for(int i = 0; i < Z.getRowCount(); i++) {
            for(int j = 0; j < Z.getColumnCount(); j++) {
                H[i][j] = Z.unsafeGet(i, j);
            }
        }
        for(int i = n - 1; i >= 0; i--) {
            H = calculateReflector(getColumn(U, i), H);
        }
        return new ComplexMatrixDense(H);
    }

    /*
     % Apply Householder transposed reflections.
    % Z = house_apply(U,X), with U from house_qr
    % computes Q'*X without actually computing Q'.
     */
    public static ComplexMatrixDense houseApplyTranspose(ComplexMatrixDense U, ComplexMatrixDense X) {
        ComplexMatrixDense Z = X.copy();
        int n = U.getColumnCount();
        Complex[][] H = new Complex[X.getRowCount()][X.getColumnCount()];
        // TODO temporary 2D array method
        for(int i = 0; i < Z.getRowCount(); i++) {
            for(int j = 0; j < Z.getColumnCount(); j++) {
                H[i][j] = Z.unsafeGet(i, j);
            }
        }
        for(int i = 0; i < n; i++) {
            H = calculateReflector(getColumn(U, i), H);
        }
        return new ComplexMatrixDense(H);
    }

    public static ComplexMatrixDense backSubstitutionSolve(ComplexMatrixDense R, ComplexMatrixDense B) {
        int n = R.getRowCount();
        int m = B.getColumnCount();

        Complex[][] X = zeros(n, m);
        for(int col = 0; col < m; col++) {
            for(int i = n - 1; i >= 0; i--) {
                X[i][col] = B.get(i, col);
                for(int j = i + 1; j < n; j++) {
                    X[i][col].subtractEquals(R.get(i, j).multiply(X[j][col]));
                }
                X[i][col].divideEquals(R.get(i, i));
            }
        }
        return new ComplexMatrixDense(X);
    }
    public static ComplexMatrixDense solve(ComplexMatrixDense A, ComplexMatrixDense B) {
        ComplexMatrixDense X = new ComplexMatrixDense(zeros(B.getRowCount(), B.getColumnCount()));

        Tuples.Tuple2<ComplexMatrixDense, ComplexMatrixDense> UR = houseQR(A);
        ComplexMatrixDense U = UR.getItem1();
        ComplexMatrixDense R = UR.getItem2();

        // Compute Y = transpose(Q) * B
        ComplexMatrixDense Y = houseApplyTranspose(U, B);

        // Solve R * X = Y;
        // Back Substitution


//            if (B.getRowCount() != rows) {
//                throw new IllegalArgumentException("Matrix row dimensions must agree.");
//            }
////            if (!this.isFullRank()) {
////                throw new RuntimeException("Matrix is rank deficient.");
////            }
//
        return backSubstitutionSolve(R, Y);
    }

    /*
    public static double[] backSubstitution(double[][] U, double[] b) {
        int n = b.length;
        double[] x = new double[n];

        for (int i = n - 1; i >= 0; i--) {
            x[i] = b[i];
            for (int j = i + 1; j < n; j++) {
                x[i] -= U[i][j] * x[j];
            }
            x[i] /= U[i][i];
        }
        return x;
    }
     */

    public ComplexQRDecompositionDense(ComplexMatrixDense matrix) {
        super(matrix);
        Complex[] data = matrix.getArrayCopy();
        _rdiag = new Complex[cols];

        for (int k = 0; k < cols; ++k) {
            double nrm = 0.0;
            // Compute 2-norm of k-th column without under/overflow.
            for (int i = k; i < rows; ++i) {
                nrm = MathETK.hypot(nrm, data[i * cols + k].abs());
            }

            if (nrm != 0.0) {
                // Form k-th Householder vector.
                if (data[k * cols + k].abs() < 0) {
                    nrm = -nrm;
                }
                for (int i = k; i < rows; i++) {
                    data[i * cols + k].divideEquals(nrm);
                }
                data[k * cols + k].addEquals(1.0);

                // Apply transformation to remaining columns.
                for (int j = k + 1; j < cols; j++) {
                    Complex s = new Complex();
                    for (int i = k; i < rows; i++) {
                        s.addEquals(data[i * cols + k].multiply(data[i * cols + j]));
                    }
                    s = s.divide(data[k * cols + k]);
                    for (int i = k; i < rows; i++) {
                        data[i * cols + j].addEquals(s.multiply(data[i * cols + k]));
                    }
                }
            }
            _rdiag[k] = Complex.fromReal(-nrm);
        }
        _data = data;
    }

    /*
     * ------------------------ Public Methods ------------------------
     */

    /**
     * Is the matrix full rank?
     *
     * @return true if R, and hence A, has full rank.
     */

//        public boolean isFullRank() {
//            for (int j = 0; j < cols; j++) {
//                if (_rdiag[j] == 0.0)
//                    return false;
//            }
//            return true;
//        }
//
//        /**
//         * Return the Householder vectors
//         *
//         * @return Lower trapezoidal matrix whose columns define the reflections
//         */
//
        public ComplexMatrixDense getH() {
            ComplexMatrixDense X = new ComplexMatrixDense(rows, cols);
            Complex[] H = X.getArray();
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    if (i >= j) {
                        H[i * cols + j] = _data[i * cols + j];
                    } else {
                        H[i * cols + j] = new Complex();
                    }
                }
            }
            return X;
        }
//
        /**
         * Return the upper triangular factor
         *
         * @return R
         */

        public ComplexMatrixDense getR() {
            ComplexMatrixDense X = new ComplexMatrixDense(cols, cols);
            Complex[] R = X.getArray();
            for (int i = 0; i < cols; i++) {
                for (int j = 0; j < cols; j++) {
                    if (i < j) {
                        R[i * cols + j] = _data[i * cols + j];
                    } else if (i == j) {
                        R[i * cols + j] = _rdiag[i];
                    } else {
                        R[i * cols + j] = new Complex();
                    }
                }
            }
            return X;
        }
//
//        /**
//         * Return the upper triangular factor
//         *
//         * @return R
//         */
//
//        public MatrixDense getRT() {
//            MatrixDense X = new MatrixDense(cols, cols);
//            double[] R = X.getArray();
//            for (int j = 0; j < cols; ++j) {
//                for (int i = 0; i < cols; ++i) {
//                    if (i > j) {
//                        R[i * cols + j] = _data[j * cols + i];
//                    } else if (i == j) {
//                        R[i * cols + j] = _rdiag[i];
//                    } else {
//                        R[i * cols + j] = 0.0;
//                    }
//                }
//            }
//            return X;
//        }
//
//        /**
//         * Generate and return the (economy-sized) orthogonal factor
//         *
//         * @return Q
//         */
//
        public ComplexMatrixDense getQThin() {
            ComplexMatrixDense X = new ComplexMatrixDense(rows, cols);
            Complex[] Q = X.getArray();
            for (int k = cols - 1; k >= 0; k--) {
                for (int i = 0; i < rows; i++) {
                    Q[i * cols + k] = new Complex();
                }
                Q[k * cols + k] = Complex.fromReal(1);
                for (int j = k; j < cols; j++) {
                    if (_data[k * cols + k].abs() != 0) {
                        Complex s = new Complex();
                        for (int i = k; i < rows; i++) {
                            s.addEquals(_data[i * cols + k].multiply(Q[i * cols + j]));
                        }
                        s = s.conj().divide(_data[k * cols + k]);
                        for (int i = k; i < rows; i++) {
                            Q[i * cols + j].addEquals(s.multiply(_data[i * cols + k]));
                        }
                    }
                }
            }
            return X;
        }


    /**
     * Generate and return the unitary orthogonal factor
     *
     * @return Q
     */
    public ComplexMatrixDense getQ() {
        // Compute Q = Q * I
        ComplexMatrixDense Q = ComplexMatrixDense.Factory.identity(rows);
        Complex[] X = Q.getArray();
        int mr = Math.min(rows, cols);
        for (int k = mr - 1; k >= 0; --k) {
            for (int j = rows - 1; j >= 0; --j) {
                Complex s = new Complex();
                for (int i = k; i < rows; i++) {
                    s.addEquals(_data[i * cols + k].multiply(X[i * rows + j]));
                }
                s = s.conj().divide(_data[k * cols + k]);
                for (int i = k; i < rows; i++) {
                    X[i * rows + j].addEquals(s.multiply(_data[i * cols + k]));
                }
            }
        }
        return Q;
    }
//
//        public MatrixDense QmultiplyX(MatrixDense X) {
//            int nx = X.getColumnCount();
//            double[] Q = X.getArray();
//            int mr = Math.min(rows, cols);
//            for (int k = mr - 1; k >= 0; --k) {
//                for (int j = nx - 1; j >= 0; --j) {
//                    double s = 0.0;
//                    for (int i = k; i < rows; i++) {
//                        s += _data[i * cols + k] * Q[i * nx + j];
//                    }
//                    s = -s / _data[k * cols + k];
//                    for (int i = k; i < rows; i++) {
//                        Q[i * nx + j] += s * _data[i * cols + k];
//                    }
//                }
//            }
//            return new MatrixDense(Q, X.getRowCount(), X.getColumnCount());
//        }
//
//        /**
//         * Generate and return the transpose of the orthogonal factor
//         *
//         * @return transpose(Q)
//         */
//        public MatrixDense getQT() {
//            // Compute Q = Q * I
//            MatrixDense Q = MatrixDense.Factory.identity(rows);
//            double[] X = Q.getArray();
//            int mr = Math.min(rows, cols);
//            for (int k = 0; k < mr; ++k) {
//                for (int j = 0; j < rows; ++j) {
//                    double s = 0.0;
//                    for (int i = k; i < rows; i++) {
//                        s += _data[i * cols + k] * X[i * rows + j];
//                    }
//                    s = -s / _data[k * cols + k];
//                    for (int i = k; i < rows; i++) {
//                        X[i * rows + j] += s * _data[i * cols + k];
//                    }
//                }
//            }
//            return Q;
//        }
//
    /**
     * Least squares solution of A*X = B
     *
     * @param B
     *            A Matrix with as many rows as A and any number of columns.
     * @return X that minimizes the two norm of Q*R*X-B.
     * @exception IllegalArgumentException
     *                Matrix row dimensions must agree.
     * @exception RuntimeException
     *                Matrix is rank deficient.
     */

//        public ComplexMatrix solve(ComplexMatrix B) {
//            if (B.getRowCount() != rows) {
//                throw new IllegalArgumentException("Matrix row dimensions must agree.");
//            }
////            if (!this.isFullRank()) {
////                throw new RuntimeException("Matrix is rank deficient.");
////            }
//
//            // Copy right hand side
//            int nx = B.getColumnCount();
//            double[] X = B.getArrayCopy();
//
//            // Compute Y = transpose(Q)*B
//            for (int k = 0; k < cols; k++) {
//                for (int j = 0; j < nx; j++) {
//                    double s = 0.0;
//                    for (int i = k; i < rows; i++) {
//                        s += _data[i * cols + k] * X[i * nx + j];
//                    }
//                    s = -s / _data[k * cols + k];
//                    for (int i = k; i < rows; i++) {
//                        X[i * nx + j] += s * _data[i * cols + k];
//                    }
//                }
//            }
//            // Solve R*X = Y;
//            for (int k = cols - 1; k >= 0; k--) {
//                for (int j = 0; j < nx; j++) {
//                    X[k * nx + j] /= _rdiag[k];
//                }
//                for (int i = 0; i < k; i++) {
//                    for (int j = 0; j < nx; j++) {
//                        X[i * nx + j] -= X[k * nx + j] * _data[i * cols + k];
//                    }
//                }
//            }
//            return (new MatrixDense(X, cols, nx).subMatrix(0, cols - 1, 0, nx - 1));
//        }
//
////	public Matrix solveTranspose(Matrix B) {
//////		if (B.getRowCount() != rows) {
//////			throw new IllegalArgumentException("Matrix row dimensions must agree.");
//////		}
////		if (!this.isFullRank()) {
////			throw new RuntimeException("Matrix is rank deficient.");
////		}
////
////		// Copy right hand side
////		int nx = B.getColumnCount();
////		double[] X = B.getArrayCopy();
////
////		// Solve RT*X = Y;
////		for (int k = cols - 1; k >= 0; k--) {
////			for (int j = 0; j < nx; j++) {
////				X[k * nx + j] /= _rdiag[k];
////			}
////			for (int i = 0; i < k; i++) {
////				for (int j = 0; j < nx; j++) {
////					X[i * nx + j] -= X[k * nx + j] * _data[i * cols + k];
////				}
////			}
////		}
////
////
////		int mr = Math.min(rows, cols);
////		for (int k = mr - 1; k >= 0; --k) {
////			for (int j = nx - 1; j >= 0; --j) {
////				double s = 0.0;
////				for (int i = k; i < rows; i++) {
////					s += _data[i * cols + k] * X[i * nx + j];
////				}
////				s = -s / _data[k * cols + k];
////				for (int i = k; i < rows; i++) {
////					X[i * nx + j] += s * _data[i * cols + k];
////				}
////			}
////		}
////		return new Matrix(X, B.getRowCount(), B.getColumnCount());
////
////		//return (new Matrix(X, cols, nx).subMatrix(0, cols - 1, 0, nx - 1));
////	}
//
//        @Override
//        public String toString() {
//            StringBuilder sb = new StringBuilder();
//            for (int i = 0; i < rows * cols; ++i) {
//                if (i > 0 && i % cols == 0) {
//                    sb.append(System.lineSeparator());
//                }
//                sb.append(String.format("%.4f", _data[i])).append(" ");
//            }
//            return sb.toString();
//        }
}
