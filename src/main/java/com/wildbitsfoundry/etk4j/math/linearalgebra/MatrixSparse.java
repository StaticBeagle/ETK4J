package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;

import java.util.Arrays;
import java.util.Iterator;

import static com.wildbitsfoundry.etk4j.math.linearalgebra.ColumnCounts.adjust;

public class MatrixSparse extends Matrix {

    /**
     * Storage for non-zero values. Only valid up to length-1.
     */
    public double[] nz_values = new double[0];
    /**
     * Length of data. Number of non-zero values in the matrix
     */
    public int nz_length;
    /**
     * Specifies which row a specific non-zero value corresponds to. If they are sorted or not with in each column
     * is specified by the {@link #indicesSorted} flag.
     */
    public int[] nz_rows = new int[0];
    /**
     * Stores the range of indexes in the non-zero lists that belong to each column. Column 'i' corresponds to
     * indexes col_idx[i] to col_idx[i+1]-1, inclusive.
     */
    public int[] col_idx;

    /**
     * Flag that's used to indicate of the row indices are sorted or not.
     */
    boolean indicesSorted = false;

//    public static double EPS = Math.pow(2.0, -52.0);

    /**
     * Constructor with a default arrayLength of zero.
     *
     * @param rows Number of rows
     * @param cols Number of columns
     */
    public MatrixSparse(int rows, int cols) {
        this(rows, cols, 0);
    }

    /**
     * Specifies shape and number of non-zero elements that can be stored.
     *
     * @param rows        Number of rows
     * @param cols        Number of columns
     * @param arrayLength Initial maximum number of non-zero elements that can be in the matrix
     */
    public MatrixSparse(int rows, int cols, int arrayLength) {
        super(rows, cols);
        if (rows < 0 || cols < 0 || arrayLength < 0)
            throw new IllegalArgumentException("Rows, columns, and arrayLength must be not be negative");
        this.nz_length = 0;
        col_idx = new int[cols + 1];
        growMaxLength(arrayLength, false);
    }

    public MatrixSparse(MatrixSparse original) {
        this(original.rows, original.cols, original.nz_length);
        nz_values = original.nz_values;
        nz_rows = original.nz_rows;
        setTo(original);
    }

    public MatrixSparse copy() {
        return new MatrixSparse(this);
    }

    public MatrixSparse createLike() {
        return new MatrixSparse(rows, cols);
    }

    public void setTo(MatrixSparse original) {
        reshape(original.rows, original.cols, original.nz_length);
        this.nz_length = original.nz_length;

        System.arraycopy(original.nz_values, 0, nz_values, 0, nz_length);
        System.arraycopy(original.nz_rows, 0, nz_rows, 0, nz_length);
        System.arraycopy(original.col_idx, 0, col_idx, 0, cols + 1);
        this.indicesSorted = original.indicesSorted;
    }

//    public void print() {
//        MatrixIO.printFancy(System.out, this, MatrixIO.DEFAULT_LENGTH);
//    }

//    public void print(String format) {
//        MatrixIO.print(System.out, this, format);
//    }

//    public void printNonZero() {
//        String format = "%d %d " + MatrixIO.DEFAULT_FLOAT_FORMAT + "\n";
//        System.out.println("Type = " + getType().name() + " , rows = " + rows + " , cols = " + cols
//                + " , nz_length = " + nz_length);
//
//        for (int col = 0; col < cols; col++) {
//            int idx0 = col_idx[col];
//            int idx1 = col_idx[col + 1];
//
//            for (int i = idx0; i < idx1; i++) {
//                int row = nz_rows[i];
//                double value = nz_values[i];
//
//                System.out.printf(format, row, col, value);
//            }
//        }
//    }

    public double[] getNonZeroValues() {
        double[] result = new double[nz_length];
        System.arraycopy(nz_values, 0, result, 0, result.length);
        return result;
    }

    public double[] getArrayDense() {
        double[] result = new double[rows * cols];
        for (int col = 0; col < cols; col++) {
            int start = col_idx[col];
            int end = col_idx[col + 1];
            for (int idx = start; idx < end; idx++) {
                int row = nz_rows[idx];
                double value = nz_values[idx];
                result[row * cols + col] = value;
            }
        }
        return result;
    }


    public boolean isAssigned(int row, int col) {
        return nz_index(row, col) >= 0;
    }

    public double get(int row, int col) {
        if (row < 0 || row >= rows || col < 0 || col >= cols)
            throw new IllegalArgumentException("Outside of matrix bounds");

        return unsafeGet(row, col);
    }

    @Override
    public double unsafeGet(int row, int col) {
        int index = nz_index(row, col);
        if (index >= 0)
            return nz_values[index];
        return 0;
    }

    public double get(int row, int col, double fallBackValue) {
        if (row < 0 || row >= rows || col < 0 || col >= cols)
            throw new IllegalArgumentException("Outside of matrix bounds");

        return unsafeGet(row, col, fallBackValue);
    }

    public double unsafeGet(int row, int col, double fallBackValue) {
        int index = nz_index(row, col);
        if (index >= 0)
            return nz_values[index];
        return fallBackValue;
    }

    /**
     * Returns the index in nz_rows for the element at (row,col) if it already exists in the matrix. If not then -1
     * is returned.
     *
     * @param row row coordinate
     * @param col column coordinate
     * @return nz_row index or -1 if the element does not exist
     */
    public int nz_index(int row, int col) {
        int col0 = col_idx[col];
        int col1 = col_idx[col + 1];

        if (this.indicesSorted) {
            return Arrays.binarySearch(nz_rows, col0, col1, row);
        } else {
            for (int i = col0; i < col1; i++) {
                if (nz_rows[i] == row) {
                    return i;
                }
            }
            return -1;
        }
    }

    public void set(int row, int col, double val) {
        if (row < 0 || row >= rows || col < 0 || col >= cols)
            throw new IllegalArgumentException("Outside of matrix bounds");

        unsafeSet(row, col, val);
    }

    /**
     * Solve system of linear equations. Three different algorithms are used depending on the shape of the matrix:
     * <pre>
     *     LU Decomposition if the matrix is squared.
     *     QR if the matrix is thin in other words it has more rows than columns. (Overdetermined system)
     *     Not supported if the matrix is short and wide in other words it has more columns than rows. (Under-determined system)
     * </pre>
     *
     * @param b The solution {@link Matrix}.
     * @return The solution to {@code Ax = b}
     */
    public MatrixSparse solve(MatrixSparse b) {
        if (rows == cols) { // Matrix is Squared
            return new LUDecompositionSparse(this).solve(b);
        } else if (rows > cols) { // Matrix is tall and narrow (Overdetermined system)
            return new QRDecompositionSparse(this).solve(b);
        } else { // Matrix is short and wide (Under-determined system)
            throw new UnsupportedOperationException("This operation is not supported for system with more columns than rows");
        }
    }

    /**
     * Solve system of linear equations. Three different algorithms are used depending on the shape of the matrix:
     * <pre>
     *     LU Decomposition if the matrix is squared.
     *     QR if the matrix is thin in other words it has more rows than columns. (Overdetermined system)
     *     Pseudo inverse * b if the matrix is short and wide in other words it has more columns than rows. (Under-determined system)
     * </pre>
     *
     * @param b The solution {@link Matrix}.
     * @return The solution to {@code Ax = b}
     */
    public MatrixSparse solve(double[] b) {
        double[][] matrix = new double[b.length][1];
        for(int i = 0; i < b.length; i++) {
            matrix[i][0] = b[i];
        }
        return solve(from2DArray(matrix));
    }

    @Override
    public double det() {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    @Override
    public LUDecompositionSparse LU() {
        return new LUDecompositionSparse(this);
    }

    @Override
    public QRDecompositionSparse QR() {
        return new QRDecompositionSparse(this);
    }

    @Override
    public CholeskyDecompositionSparse Chol() { return new CholeskyDecompositionSparse(this); }

    @Override
    public boolean isEmpty() {
        return rows == 0 && cols == 0;
    }

    /**
     * Inverse of the {@code Matrix}.
     *
     * @return {@code A<sup>-1</sup>}.
     */
    public MatrixSparse inv() {
        return this.solve(MatrixSparse.Factory.identity(rows));
    }

    @Override
    public void unsafeSet(int row, int col, double val) {
        int index = nz_index(row, col);
        if (index >= 0) {
            nz_values[index] = val;
        } else {

            int idx0 = col_idx[col];
            int idx1 = col_idx[col + 1];

            // determine the index the new element should be inserted at. This is done to keep it sorted if
            // it was already sorted
            for (index = idx0; index < idx1; index++) {
                if (row < nz_rows[index]) {
                    break;
                }
            }

            // shift all the col_idx after this point by 1
            for (int i = col + 1; i <= cols; i++) {
                col_idx[i]++;
            }

            // if it's already at the maximum array length grow the arrays
            if (nz_length >= nz_values.length)
                growMaxLength(nz_length * 2 + 1, true);

            // shift everything by one
            for (int i = nz_length; i > index; i--) {
                nz_rows[i] = nz_rows[i - 1];
                nz_values[i] = nz_values[i - 1];
            }
            nz_rows[index] = row;
            nz_values[index] = val;
            nz_length++;
        }
    }

    public void remove(int row, int col) {
        int index = nz_index(row, col);

        if (index < 0) // it's not in the nz structure
            return;

        // shift all the col_idx after this point by -1
        for (int i = col + 1; i <= cols; i++) {
            col_idx[i]--;
        }

        nz_length--;
        for (int i = index; i < nz_length; i++) {
            nz_rows[i] = nz_rows[i + 1];
            nz_values[i] = nz_values[i + 1];
        }
    }

    public void zero() {
        Arrays.fill(col_idx, 0, cols + 1, 0);
        nz_length = 0;
        indicesSorted = false; // see justification in reshape
    }

    public MatrixSparse create(int rows, int cols) {
        return new MatrixSparse(rows, cols);
    }

    public int getNonZeroLength() {
        return nz_length;
    }

    public void reshape(int rows, int cols, int arrayLength) {
        if (rows < 0 || cols < 0 || arrayLength < 0)
            throw new IllegalArgumentException("Rows, columns, and arrayLength must be not be negative");

        // OK so technically it is sorted, but forgetting to correctly set this flag is a common mistake so
        // decided to be conservative and mark it as unsorted so that stuff doesn't blow up
        this.indicesSorted = false;
        this.rows = rows;
        this.cols = cols;
        growMaxLength(arrayLength, false);
        this.nz_length = 0;

        if (cols + 1 > col_idx.length) {
            col_idx = new int[cols + 1];
        } else {
            Arrays.fill(col_idx, 0, cols + 1, 0);
        }
    }


    public void reshape(int rows, int cols) {
        reshape(rows, cols, 0);
    }


    public void shrinkArrays() {
        if (nz_length < nz_values.length) {
            double[] tmp_values = new double[nz_length];
            int[] tmp_rows = new int[nz_length];

            System.arraycopy(this.nz_values, 0, tmp_values, 0, nz_length);
            System.arraycopy(this.nz_rows, 0, tmp_rows, 0, nz_length);

            this.nz_values = tmp_values;
            this.nz_rows = tmp_rows;
        }
    }

    /**
     * Increases the maximum size of the data array so that it can store sparse data up to 'length'. The class
     * parameter nz_length is not modified by this function call.
     *
     * @param arrayLength   Desired maximum length of sparse data
     * @param preserveValue If true the old values will be copied into the new arrays. If false that step will be skipped.
     */
    public void growMaxLength(int arrayLength, boolean preserveValue) {
        if (arrayLength < 0)
            throw new IllegalArgumentException("Negative array length. Overflow?");

        // NOTE: The code below has been (experimentally) commented out. A situation arose where we wanted to exceed
        //       the max physical size, which would then be corrected later on.

        // see if multiplying rows*cols will cause an overflow. If it won't then pick the smaller of the two
//        if( rows != 0 && cols <= Integer.MAX_VALUE / rows ) {
//            // save the user from themselves
//            arrayLength = Math.min(rows*cols, arrayLength);
//        }
        if (arrayLength > this.nz_values.length) {
            double[] data = new double[arrayLength];
            int[] row_idx = new int[arrayLength];

            if (preserveValue) {
                System.arraycopy(this.nz_values, 0, data, 0, this.nz_length);
                System.arraycopy(this.nz_rows, 0, row_idx, 0, this.nz_length);
            }

            this.nz_values = data;
            this.nz_rows = row_idx;
        }
    }

    /**
     * Increases the maximum number of columns in the matrix.
     *
     * @param desiredColumns Desired number of columns.
     * @param preserveValue  If the array needs to be expanded should it copy the previous values?
     */
    public void growMaxColumns(int desiredColumns, boolean preserveValue) {
        if (col_idx.length < desiredColumns + 1) {
            int[] c = new int[desiredColumns + 1];
            if (preserveValue)
                System.arraycopy(col_idx, 0, c, 0, col_idx.length);
            col_idx = c;
        }
    }

    /**
     * Given the histogram of columns compute the col_idx for the matrix. nz_length is automatically set and
     * nz_values will grow if needed.
     *
     * @param histogram histogram of column values in the sparse matrix. modified, see above.
     */
    public void histogramToStructure(int[] histogram) {
        col_idx[0] = 0;
        int index = 0;
        for (int i = 1; i <= cols; i++) {
            col_idx[i] = index += histogram[i - 1];
        }
        nz_length = index;
        growMaxLength(nz_length, false);
        if (col_idx[cols] != nz_length)
            throw new RuntimeException("Egads");
    }

    /**
     * Copies the non-zero structure of orig into "this"
     *
     * @param orig Matrix who's structure is to be copied
     */
    public void copyStructure(MatrixSparse orig) {
        reshape(orig.rows, orig.cols, orig.nz_length);
        this.nz_length = orig.nz_length;
        System.arraycopy(orig.col_idx, 0, col_idx, 0, orig.cols + 1);
        System.arraycopy(orig.nz_rows, 0, nz_rows, 0, orig.nz_length);
    }

    /**
     * If the indices has been sorted or not
     *
     * @return true if sorted or false if not sorted
     */
    public boolean isIndicesSorted() {
        return indicesSorted;
    }

    /**
     * Returns true if number of non-zero elements is the maximum size
     *
     * @return true if no more non-zero elements can be added
     */
    public boolean isFull() {
        return nz_length == rows * cols;
    }

    public static MatrixSparse from2DArray(double[][] array) {
        return from2DArray(array, ConstantsETK.DOUBLE_EPS);
    }

    public static MatrixSparse from2DArray(double[][] array, double tol) {
        int rows = array.length;
        int cols = array[0].length;

        MatrixSparse matrixSparse = new MatrixSparse(rows, cols);

        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                if(Math.abs(array[i][j]) > tol) {
                    matrixSparse.unsafeSet(i, j, array[i][j]);
                }
            }
        }
        return matrixSparse;
    }

    /**
     * Value of an element in a sparse matrix
     */
    public static class CoordinateRealValue {
        /**
         * The coordinate
         */
        public int row, col;
        /**
         * The value of the coordinate
         */
        public double value;
    }

    public Iterator<CoordinateRealValue> createCoordinateIterator() {
        return new Iterator<CoordinateRealValue>() {
            final CoordinateRealValue coordinate = new CoordinateRealValue();
            int nz_index = 0; // the index of the non-zero value and row
            int column = 0; // which column it's in

            {
                incrementColumn();
            }


            public boolean hasNext() {
                return nz_index < nz_length;
            }


            public CoordinateRealValue next() {
                coordinate.row = nz_rows[nz_index];
                coordinate.col = column;
                coordinate.value = nz_values[nz_index];
                nz_index++;
                incrementColumn();
                return coordinate;
            }

            private void incrementColumn() {
                while (column + 1 <= cols && nz_index >= col_idx[column + 1]) {
                    column++;
                }
            }
        };
    }

    /**
     * Performs matrix addition:<br>
     * C = this + B
     *
     * @param B Matrix
     * @return this + B
     */
    public MatrixSparse add(MatrixSparse B) {
        MatrixSparse outputC = new MatrixSparse(rows, cols);
        if (rows != B.rows || cols != B.cols)
            throw new IllegalArgumentException("Inconsistent matrix shapes.");
        outputC = reshapeOrDeclare(outputC, this, rows, cols);

        addOp(this, B, outputC);

        return outputC;
    }

    public static MatrixSparse reshapeOrDeclare(MatrixSparse target, MatrixSparse reference, int rows, int cols) {
        if (target == null)
            return reference.create(rows, cols);
        else if (target.rows != rows || target.cols != cols)
            target.reshape(rows, cols);
        return target;
    }

    /**
     * Performs matrix addition:<br>
     * C = A + B
     *
     * @param A Matrix
     * @param B Matrix
     * @param C Output matrix.
     */
    private static void addOp(MatrixSparse A, MatrixSparse B, MatrixSparse C) {
        double[] x = adjust(new DGrowArray(), A.rows);
        int[] w = adjust(new IGrowArray(), A.rows, A.rows);

        C.indicesSorted = false;
        C.nz_length = 0;

        for (int col = 0; col < A.cols; col++) {
            C.col_idx[col] = C.nz_length;

            multiplyAddColA(A, col, 1, C, col + 1, x, w);
            multiplyAddColA(B, col, 1, C, col + 1, x, w);

            // take the values in the dense vector 'x' and put them into 'C'
            int idxC0 = C.col_idx[col];
            int idxC1 = C.col_idx[col + 1];

            for (int i = idxC0; i < idxC1; i++) {
                C.nz_values[i] = x[C.nz_rows[i]];
            }
        }
        C.col_idx[A.cols] = C.nz_length;
    }

    /**
     * Performs matrix addition:<br>
     * C = this - B
     *
     * @param B Matrix
     * @return this - B
     */
    public MatrixSparse subtract(MatrixSparse B) {
        MatrixSparse outputC = new MatrixSparse(rows, cols);
        if (rows != B.rows || cols != B.cols)
            throw new IllegalArgumentException("Inconsistent matrix shapes.");
        outputC = reshapeOrDeclare(outputC, this, rows, cols);

        subtractOp(this, B, outputC);

        return outputC;
    }

    /**
     * Performs A - B
     * @param A Matrix
     * @param B Matrix
     * @param C Output matrix.
     */
    private static void subtractOp(MatrixSparse A, MatrixSparse B, MatrixSparse C) {
        double[] x = adjust(new DGrowArray(), A.rows);
        int[] w = adjust(new IGrowArray(), A.rows, A.rows);

        C.indicesSorted = false;
        C.nz_length = 0;

        for (int col = 0; col < A.cols; col++) {
            C.col_idx[col] = C.nz_length;

            multiplyAddColA(A, col, 1, C, col + 1, x, w);
            multiplyAddColA(B, col, -1, C, col + 1, x, w);

            // take the values in the dense vector 'x' and put them into 'C'
            int idxC0 = C.col_idx[col];
            int idxC1 = C.col_idx[col + 1];

            for (int i = idxC0; i < idxC1; i++) {
                C.nz_values[i] = x[C.nz_rows[i]];
            }
        }
        C.col_idx[A.cols] = C.nz_length;
    }

    /**
     * Performs matrix multiplication. C = this * B
     *
     * @param B Matrix
     * @return this * B
     */
    public MatrixSparse multiply(MatrixSparse B) {
        MatrixSparse C = new MatrixSparse(this.rows, B.cols);
        multiplyOp(this, B, C);
        return C;
    }

    /**
     * Performs matrix multiplication. C = this * scalar
     *
     * @param scalar value
     * @return this * scalar
     */
    public MatrixSparse multiply(double scalar) {
        MatrixSparse C = new MatrixSparse(rows, cols);
        for (int col = 0; col < cols; col++) {
            int start = col_idx[col];
            int end = col_idx[col + 1];
            for (int idx = start; idx < end; idx++) {
                int row = nz_rows[idx];
                double value = nz_values[idx];
                C.unsafeSet(row, col, value * scalar);
            }
        }
        return C;
    }

    /**
     * Performs matrix multiplication. C = A * B
     *
     * @param A Matrix
     * @param B Matrix
     * @param C Storage for results. Array size is increased if needed.
     */
    private static void multiplyOp(MatrixSparse A, MatrixSparse B, MatrixSparse C) {

        double[] x = adjust(new DGrowArray(), A.rows);
        int[] w = adjust(new IGrowArray(), A.rows, A.rows);

        C.growMaxLength(A.nz_length + B.nz_length, false);
        C.indicesSorted = false;
        C.nz_length = 0;

        // C(i,j) = sum_k A(i,k) * B(k,j)
        int idx0 = B.col_idx[0];
        for (int bj = 1; bj <= B.cols; bj++) {
            int colB = bj - 1;
            int idx1 = B.col_idx[bj];
            C.col_idx[bj] = C.nz_length;

            if (idx0 == idx1) {
                continue;
            }

            // C(:,j) = sum_k A(:,k)*B(k,j)
            for (int bi = idx0; bi < idx1; bi++) {
                int rowB = B.nz_rows[bi];
                double valB = B.nz_values[bi];  // B(k,j)  k=rowB j=colB

                multiplyAddColA(A, rowB, valB, C, colB + 1, x, w);
            }

            // take the values in the dense vector 'x' and put them into 'C'
            int idxC0 = C.col_idx[colB];
            int idxC1 = C.col_idx[colB + 1];

            for (int i = idxC0; i < idxC1; i++) {
                C.nz_values[i] = x[C.nz_rows[i]];
            }

            idx0 = idx1;
        }
    }

    /*
     * Performs the operation x = x + A(:,i)*alpha
     *
     * <p>NOTE: This is the same as cs_scatter() in csparse.</p>
     */
    private static void multiplyAddColA(MatrixSparse A, int colA,
                                        double alpha,
                                        MatrixSparse C, int mark,
                                        double[] x, int[] w) {
        int idxA0 = A.col_idx[colA];
        int idxA1 = A.col_idx[colA + 1];

        for (int j = idxA0; j < idxA1; j++) {
            int row = A.nz_rows[j];

            if (w[row] < mark) {
                if (C.nz_length >= C.nz_rows.length) {
                    C.growMaxLength(C.nz_length * 2 + 1, true);
                }

                w[row] = mark;
                C.nz_rows[C.nz_length] = row;
                C.col_idx[mark] = ++C.nz_length;
                x[row] = A.nz_values[j] * alpha;
            } else {
                x[row] += A.nz_values[j] * alpha;
            }
        }
    }

    /**
     * Perform matrix transpose
     *
     * @return The transposed matrix
     */
    public MatrixSparse transpose() {
        MatrixSparse A_t = reshapeOrDeclare(null, cols, rows, nz_length);
        transposeOp(this, A_t);
        return A_t;
    }

    public static MatrixSparse reshapeOrDeclare(MatrixSparse target, int rows, int cols, int nz_length) {
        if (target == null)
            return new MatrixSparse(rows, cols, nz_length);
        else
            target.reshape(rows, cols, nz_length);
        return target;
    }

    private static void transposeOp(MatrixSparse A, MatrixSparse C) {
        int[] work = adjust((IGrowArray) null, A.cols, A.rows);
        C.reshape(A.cols, A.rows, A.nz_length);

        // compute the histogram for each row in 'a'
        for (int j = 0; j < A.nz_length; j++) {
            work[A.nz_rows[j]]++;
        }

        // construct col_idx in the transposed matrix
        C.histogramToStructure(work);
        System.arraycopy(C.col_idx, 0, work, 0, C.cols);

        // fill in the row indexes
        int idx0 = A.col_idx[0];
        for (int j = 1; j <= A.cols; j++) {
            final int col = j - 1;
            final int idx1 = A.col_idx[j];
            for (int i = idx0; i < idx1; i++) {
                int row = A.nz_rows[i];
                int index = work[row]++;
                C.nz_rows[index] = col;
                C.nz_values[index] = A.nz_values[i];
            }
            idx0 = idx1;
        }
    }

    public MatrixDense toDense() {
        MatrixDense matrixDense = new MatrixDense(rows, cols);
        for (int col = 0; col < cols; col++) {
            int start = col_idx[col];
            int end = col_idx[col + 1];
            for (int idx = start; idx < end; idx++) {
                int row = nz_rows[idx];
                double value = nz_values[idx];
                matrixDense.unsafeSet(row, col, value);
            }
        }
        return matrixDense;
    }

    @Override
    public String toString() {
        return toDense().toString()
                .replaceAll(" {6}0.000 {6}", "                 ")
                .replaceAll("0.000 {6}", "           ");
    }

    public static class Factory {
        private Factory() {}

        /**
         * Identity {@code Matrix}.
         * @param rows The number of rows.
         * @param cols The number of columns.
         * @return {@code identity(rows, cols)}.
         */
        public static MatrixSparse identity(int rows, int cols) {
            MatrixSparse identity = new MatrixSparse(rows, cols);
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    if (i == j) {
                        identity.unsafeSet(i, j, 1);
                    }
                }
            }
            return identity;
        }

        /**
         * Identity {@code Matrix.}
         * @param n The number of rows and columns.
         * @return {@code identity(n, n)}.
         */
        public static MatrixSparse identity(int n) {
            return Factory.identity(n, n);
        }
    }
}
