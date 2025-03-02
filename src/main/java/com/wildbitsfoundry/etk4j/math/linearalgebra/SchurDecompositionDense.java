package com.wildbitsfoundry.etk4j.math.linearalgebra;

/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import com.wildbitsfoundry.etk4j.constant.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.signal.filters.MaximumNumberOfIterationsExceededException;


/**
 * Class transforming a general real matrix to Schur form.
 * <p>A m × m matrix A can be written as the product of three matrices: A = P
 * × T × P<sup>T</sup> with P an orthogonal matrix and T an quasi-triangular
 * matrix. Both P and T are m × m matrices.</p>
 * <p>Transformation to Schur form is often not a goal by itself, but it is an
 * intermediate step in more general decomposition algorithms like
 * {@link EigenvalueDecompositionDense eigen decomposition}. This class is therefore
 * intended for internal use by the library and is not public. As a consequence
 * of this explicitly limited scope, many methods directly returns references to
 * internal arrays, not copies.</p>
 * <p>This class is based on the method hqr2 in class EigenvalueDecomposition
 * from the <a href="http://math.nist.gov/javanumerics/jama/">JAMA library</a>.
 * */
public class SchurDecompositionDense {

    /**
     * Maximum allowed iterations for convergence of the transformation.
     */
    private static final int MAX_ITERATIONS = 100;

    /**
     * U matrix.
     */
    private final double[][] matrixU;
    /**
     * T matrix.
     */
    private final double[][] matrixT;
    /**
     * Cached value of P.
     */
    private MatrixDense cachedU;
    /**
     * Cached value of T.
     */
    private MatrixDense cachedT;
    /**
     * Cached value of PT.
     */
    private MatrixDense cachedUt;

    /**
     * Epsilon criteria taken from JAMA code (originally was 2^-52).
     */
    private final double epsilon = ConstantsETK.DOUBLE_EPS;

    /**
     * Build the transformation to Schur form of a general real matrix.
     *
     * @param matrix matrix to transform
     * @throws NonSquareMatrixException if the matrix is not square
     */
    SchurDecompositionDense(final MatrixDense matrix) {
        if (!matrix.isSquare()) {
            throw new NonSquareMatrixException("Matrix must be a square matrix");
        }

        HessembergDecompositionDense hess = new HessembergDecompositionDense(matrix);
        matrixT = hess.getH().getAs2DArray();
        matrixU = hess.getU().getAs2DArray();
        cachedT = null;
        cachedU = null;
        cachedUt = null;

        // transform matrix
        transform();
    }

    /**
     * Returns the matrix P of the transform.
     * <p>P is an orthogonal matrix, i.e. its inverse is also its transpose.
     *
     * @return the P matrix
     */
    public MatrixDense getU() {
        if (cachedU == null) {
            cachedU = MatrixDense.from2DArray(matrixU);
        }
        return cachedU;
    }

    /**
     * Returns the transpose of the matrix P of the transform.
     * <p>P is an orthogonal matrix, i.e. its inverse is also its transpose.
     *
     * @return the transpose of the P matrix
     */
    public MatrixDense getUT() {
        if (cachedUt == null) {
            cachedUt = getU().transpose();
        }

        // return the cached matrix
        return cachedUt;
    }

    /**
     * Returns the quasi-triangular Schur matrix T of the transform.
     *
     * @return the T matrix
     */
    public MatrixDense getT() {
        if (cachedT == null) {
            cachedT = MatrixDense.from2DArray(matrixT);
        }

        // return the cached matrix
        return cachedT;
    }

    /**
     * Transform original matrix to Schur form.
     */
    private void transform() {
        final int n = matrixT.length;

        // compute matrix norm
        final double norm = getNorm();

        // shift information
        final ShiftInfo shift = new ShiftInfo();

        // Outer loop over eigenvalue index
        int iteration = 0;
        int iu = n - 1;
        while (iu >= 0) {

            // Look for single small sub-diagonal element
            final int il = findSmallSubDiagonalElement(iu, norm);

            // Check for convergence
            if (il == iu) {
                // One root found
                matrixT[iu][iu] += shift.exShift;
                iu--;
                iteration = 0;
            } else if (il == iu - 1) {
                // Two roots found
                double p = (matrixT[iu - 1][iu - 1] - matrixT[iu][iu]) / 2.0;
                double q = p * p + matrixT[iu][iu - 1] * matrixT[iu - 1][iu];
                matrixT[iu][iu] += shift.exShift;
                matrixT[iu - 1][iu - 1] += shift.exShift;

                if (q >= 0) {
                    double z = Math.sqrt(Math.abs(q));
                    if (p >= 0) {
                        z = p + z;
                    } else {
                        z = p - z;
                    }
                    final double x = matrixT[iu][iu - 1];
                    final double s = Math.abs(x) + Math.abs(z);
                    p = x / s;
                    q = z / s;
                    final double r = Math.sqrt(p * p + q * q);
                    p /= r;
                    q /= r;

                    // Row modification
                    for (int j = iu - 1; j < n; j++) {
                        z = matrixT[iu - 1][j];
                        matrixT[iu - 1][j] = q * z + p * matrixT[iu][j];
                        matrixT[iu][j] = q * matrixT[iu][j] - p * z;
                    }

                    // Column modification
                    for (int i = 0; i <= iu; i++) {
                        z = matrixT[i][iu - 1];
                        matrixT[i][iu - 1] = q * z + p * matrixT[i][iu];
                        matrixT[i][iu] = q * matrixT[i][iu] - p * z;
                    }

                    // Accumulate transformations
                    for (int i = 0; i <= n - 1; i++) {
                        z = matrixU[i][iu - 1];
                        matrixU[i][iu - 1] = q * z + p * matrixU[i][iu];
                        matrixU[i][iu] = q * matrixU[i][iu] - p * z;
                    }
                }
                iu -= 2;
                iteration = 0;
            } else {
                // No convergence yet
                computeShift(il, iu, iteration, shift);

                // stop transformation after too many iterations
                if (++iteration > MAX_ITERATIONS) {
                    throw new MaximumNumberOfIterationsExceededException("Maximum number of iterations exceeded");
                }

                // the initial houseHolder vector for the QR step
                final double[] hVec = new double[3];

                final int im = initQRStep(il, iu, shift, hVec);
                performDoubleQRStep(il, im, iu, shift, hVec);
            }
        }
    }

    /**
     * Computes the L1 norm of the (quasi-)triangular matrix T.
     *
     * @return the L1 norm of matrix T
     */
    private double getNorm() {
        double norm = 0.0;
        for (int i = 0; i < matrixT.length; i++) {
            // as matrix T is (quasi-)triangular, also take the sub-diagonal element into account
            for (int j = Math.max(i - 1, 0); j < matrixT.length; j++) {
                norm += Math.abs(matrixT[i][j]);
            }
        }
        return norm;
    }

    /**
     * Find the first small sub-diagonal element and returns its index.
     *
     * @param startIdx the starting index for the search
     * @param norm     the L1 norm of the matrix
     * @return the index of the first small sub-diagonal element
     */
    private int findSmallSubDiagonalElement(final int startIdx, final double norm) {
        int l = startIdx;
        while (l > 0) {
            double s = Math.abs(matrixT[l - 1][l - 1]) + Math.abs(matrixT[l][l]);
            if (s == 0.0) {
                s = norm;
            }
            if (Math.abs(matrixT[l][l - 1]) < epsilon * s) {
                break;
            }
            l--;
        }
        return l;
    }

    /**
     * Compute the shift for the current iteration.
     *
     * @param l         the index of the small sub-diagonal element
     * @param idx       the current eigenvalue index
     * @param iteration the current iteration
     * @param shift     holder for shift information
     */
    private void computeShift(final int l, final int idx, final int iteration, final ShiftInfo shift) {
        // Form shift
        shift.x = matrixT[idx][idx];
        shift.y = shift.w = 0.0;
        if (l < idx) {
            shift.y = matrixT[idx - 1][idx - 1];
            shift.w = matrixT[idx][idx - 1] * matrixT[idx - 1][idx];
        }

        // Wilkinson's original ad hoc shift
        if (iteration == 10) {
            shift.exShift += shift.x;
            for (int i = 0; i <= idx; i++) {
                matrixT[i][i] -= shift.x;
            }
            final double s = Math.abs(matrixT[idx][idx - 1]) + Math.abs(matrixT[idx - 1][idx - 2]);
            shift.x = 0.75 * s;
            shift.y = 0.75 * s;
            shift.w = -0.4375 * s * s;
        }

        // MATLAB's new ad hoc shift
        if (iteration == 30) {
            double s = (shift.y - shift.x) / 2.0;
            s = s * s + shift.w;
            if (s > 0.0) {
                s = Math.sqrt(s);
                if (shift.y < shift.x) {
                    s = -s;
                }
                s = shift.x - shift.w / ((shift.y - shift.x) / 2.0 + s);
                for (int i = 0; i <= idx; i++) {
                    matrixT[i][i] -= s;
                }
                shift.exShift += s;
                shift.x = shift.y = shift.w = 0.964;
            }
        }
    }

    /**
     * Initialize the householder vectors for the QR step.
     *
     * @param il    the index of the small sub-diagonal element
     * @param iu    the current eigenvalue index
     * @param shift shift information holder
     * @param hVec  the initial houseHolder vector
     * @return the start index for the QR step
     */
    private int initQRStep(int il, final int iu, final ShiftInfo shift, double[] hVec) {
        // Look for two consecutive small sub-diagonal elements
        int im = iu - 2;
        while (im >= il) {
            final double z = matrixT[im][im];
            final double r = shift.x - z;
            double s = shift.y - z;
            hVec[0] = (r * s - shift.w) / matrixT[im + 1][im] + matrixT[im][im + 1];
            hVec[1] = matrixT[im + 1][im + 1] - z - r - s;
            hVec[2] = matrixT[im + 2][im + 1];

            if (im == il) {
                break;
            }

            final double lhs = Math.abs(matrixT[im][im - 1]) * (Math.abs(hVec[1]) + Math.abs(hVec[2]));
            final double rhs = Math.abs(hVec[0]) * (Math.abs(matrixT[im - 1][im - 1]) +
                    Math.abs(z) +
                    Math.abs(matrixT[im + 1][im + 1]));

            if (lhs < epsilon * rhs) {
                break;
            }
            im--;
        }

        return im;
    }

    /**
     * Perform a double QR step involving rows l:idx and columns m:n
     *
     * @param il    the index of the small sub-diagonal element
     * @param im    the start index for the QR step
     * @param iu    the current eigenvalue index
     * @param shift shift information holder
     * @param hVec  the initial houseHolder vector
     */
    private void performDoubleQRStep(final int il, final int im, final int iu,
                                     final ShiftInfo shift, final double[] hVec) {

        final int n = matrixT.length;
        double p = hVec[0];
        double q = hVec[1];
        double r = hVec[2];

        for (int k = im; k <= iu - 1; k++) {
            boolean notlast = k != (iu - 1);
            if (k != im) {
                p = matrixT[k][k - 1];
                q = matrixT[k + 1][k - 1];
                r = notlast ? matrixT[k + 2][k - 1] : 0.0;
                shift.x = Math.abs(p) + Math.abs(q) + Math.abs(r);
                if (MathETK.isClose(shift.x, 0.0, epsilon, 0)) {
                    continue;
                }
                p /= shift.x;
                q /= shift.x;
                r /= shift.x;
            }
            double s = Math.sqrt(p * p + q * q + r * r);
            if (p < 0.0) {
                s = -s;
            }
            if (s != 0.0) {
                if (k != im) {
                    matrixT[k][k - 1] = -s * shift.x;
                } else if (il != im) {
                    matrixT[k][k - 1] = -matrixT[k][k - 1];
                }
                p += s;
                shift.x = p / s;
                shift.y = q / s;
                double z = r / s;
                q /= p;
                r /= p;

                // Row modification
                for (int j = k; j < n; j++) {
                    p = matrixT[k][j] + q * matrixT[k + 1][j];
                    if (notlast) {
                        p += r * matrixT[k + 2][j];
                        matrixT[k + 2][j] -= p * z;
                    }
                    matrixT[k][j] -= p * shift.x;
                    matrixT[k + 1][j] -= p * shift.y;
                }

                // Column modification
                for (int i = 0; i <= Math.min(iu, k + 3); i++) {
                    p = shift.x * matrixT[i][k] + shift.y * matrixT[i][k + 1];
                    if (notlast) {
                        p += z * matrixT[i][k + 2];
                        matrixT[i][k + 2] -= p * r;
                    }
                    matrixT[i][k] -= p;
                    matrixT[i][k + 1] -= p * q;
                }

                // Accumulate transformations
                final int high = matrixT.length - 1;
                for (int i = 0; i <= high; i++) {
                    p = shift.x * matrixU[i][k] + shift.y * matrixU[i][k + 1];
                    if (notlast) {
                        p += z * matrixU[i][k + 2];
                        matrixU[i][k + 2] -= p * r;
                    }
                    matrixU[i][k] -= p;
                    matrixU[i][k + 1] -= p * q;
                }
            }  // (s != 0)
        }  // k loop

        // clean up pollution due to round-off errors
        for (int i = im + 2; i <= iu; i++) {
            matrixT[i][i - 2] = 0.0;
            if (i > im + 2) {
                matrixT[i][i - 3] = 0.0;
            }
        }
    }

    /**
     * Internal data structure holding the current shift information.
     * Contains variable names as present in the original JAMA code.
     */
    private static class ShiftInfo {
        // CHECKSTYLE: stop all

        /**
         * x shift info
         */
        double x;
        /**
         * y shift info
         */
        double y;
        /**
         * w shift info
         */
        double w;
        /**
         * Indicates an exceptional shift.
         */
        double exShift;

        // CHECKSTYLE: resume all
    }
}
