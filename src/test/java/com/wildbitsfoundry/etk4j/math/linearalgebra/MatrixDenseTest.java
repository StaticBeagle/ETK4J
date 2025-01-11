package com.wildbitsfoundry.etk4j.math.linearalgebra;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;

public class MatrixDenseTest {

    @Test
    public void testMagic() {

        MatrixDense sol = new MatrixDense(1, 1, 1);
        MatrixDense magic = MatrixDense.Factory.magic(1);
        assertEquals(sol, magic);

        sol = MatrixDense.Factory.empty();
        magic = MatrixDense.Factory.magic(2);
        assertEquals(sol, magic);

        double[][] values3by3 = {{8.0, 1.0, 6.0}, {3.0, 5.0, 7.0}, {4.0, 9.0, 2.0}};
        sol = new MatrixDense(values3by3);
        magic = MatrixDense.Factory.magic(3);
        assertEquals(sol, magic);

        double[][] values4by4 = {{16.0, 2.0, 3.0, 13.0}, {5.0, 11.0, 10.0, 8.0}, {9.0, 7.0, 6.0, 12.0},
                {4.0, 14.0, 15.0, 1.0}};
        sol = new MatrixDense(values4by4);
        magic = MatrixDense.Factory.magic(4);
        assertEquals(sol, magic);

        double[][] values5by5 = {{17.0, 24.0, 1.0, 8.0, 15.0}, {23.0, 5.0, 7.0, 14.0, 16.0},
                {4.0, 6.0, 13.0, 20.0, 22.0}, {10.0, 12.0, 19.0, 21.0, 3.0}, {11.0, 18.0, 25.0, 2.0, 9.0}};
        sol = new MatrixDense(values5by5);
        magic = MatrixDense.Factory.magic(5);
        assertEquals(sol, magic);

        double[][] values6by6 = {{35.0, 1.0, 6.0, 26.0, 19.0, 24.0}, {3.0, 32.0, 7.0, 21.0, 23.0, 25.0},
                {31.0, 9.0, 2.0, 22.0, 27.0, 20.0}, {8.0, 28.0, 33.0, 17.0, 10.0, 15.0},
                {30.0, 5.0, 34.0, 12.0, 14.0, 16.0}, {4.0, 36.0, 29.0, 13.0, 18.0, 11.0}};
        sol = new MatrixDense(values6by6);
        magic = MatrixDense.Factory.magic(6);
        assertEquals(sol, magic);
    }

    @Test
    public void testPinv() {
        double[][] values = {{1, 2}, {3, 4}, {5, 6}};
        MatrixDense matrix = new MatrixDense(values);

        double[] solution = {-1.3333333333333317, -0.3333333333333325, 0.6666666666666661, 1.0833333333333321,
                0.33333333333333265, -0.4166666666666662};
        assertEquals(new MatrixDense(solution, 2, 3), matrix.pinv());
    }

    @Test
    public void testExpm() {
        double[][] data = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
        MatrixDense matrix = new MatrixDense(data);

        double[] solution = {1118906.6994131864, 1374815.0629358063, 1630724.4264584277, 2533881.0418989714,
                3113415.0313805523, 3692947.020862136, 3948856.3843847574, 4852012.9998253, 5755170.615265846};
        assertEquals(new MatrixDense(solution, 3, 3), matrix.expm());
    }

    @Test
    public void testSetRow() {
        MatrixDense m = MatrixDense.Factory.magic(3);
        double[] row = {-1, -1, -1};
        m.setRow(1, row);
        assertArrayEquals(row, m.getRow(1), 1e-12);
    }

    @Test
    public void testGetCol() {
        MatrixDense m = MatrixDense.Factory.magic(3);
        double[] row = {8, 3, 4};
        assertArrayEquals(row, m.getCol(0), 1e-12);
    }

    @Test
    public void testHessembergDecomposition() {
        double[][] matrix = {
                {65, 35, 40, 69},
                {99, 64, 37, 2},
                {39, 48, 35, 90},
                {30, 93, 87, 17}
        };

        MatrixDense A = MatrixDense.from2DArray(matrix);
        HessembergDecompositionDense hess = new HessembergDecompositionDense(A);

        double[] expectedH = {65.0, -64.17727312578465, -58.72784396012248, 4.279948356457808, -110.55315463612966,
                123.81148748159052, 17.486023979786594, 22.14752145860223, 0.0, 100.60648832932314, 36.36015423097008,
                26.39970581394131, 0.0, 0.0, 58.20341156591577, -44.17164171256066};
        assertArrayEquals(expectedH, hess.getH().getArray(), 1e-12);
        double[] expectedU = {1.0, 0.0, 0.0, 0.0, 0.0, -0.8954968343132741, 0.39724801290719464, 0.2006973741138373,
                0.0, -0.3527714801840171, -0.3585884935028961, -0.8642722806477718, 0.0, -0.2713626770646286,
                -0.8447534010399773, 0.46125263026644026};
        assertArrayEquals(expectedU, hess.getU().getArray(), 1e-12);
    }

    @Test
    public void testSchurDecomposition() {
        double[][] matrix = {
                {65, 35, 40, 69},
                {99, 64, 37, 2},
                {39, 48, 35, 90},
                {30, 93, 87, 17}
        };

        MatrixDense A = MatrixDense.from2DArray(matrix);
        SchurDecompositionDense schur = new SchurDecompositionDense(A);

        double[] expectedT = {211.53821949564204, 27.963104879902804, 11.529192222479173, -11.83821931259974,
                -2.8398992587956425E-29, 27.374680489564316, 64.83795403092125, 3.6006930114899935, 0.0,
                -31.22064867637267, -13.419034896325854, 61.10465539906722, 0.0, 0.0, 2.2836324137572356E-24,
                -44.49386508888031};
        assertArrayEquals(expectedT, schur.getT().getArray(), 1e-12);
        double[] expectedP = {0.49782038512887483, 0.15442410294358475, 0.7718547804323013, 0.3640992426578396,
                0.4680103341616455, 0.7666501603471897, -0.30669802271073127, -0.3148812182758139, 0.505725691055896,
                -0.5430518840085535, 0.0954550851000274, -0.6634941623023844, 0.5266713554713984, -0.3057701413553275,
                -0.5486937647884095, 0.5727626528185856};
        assertArrayEquals(expectedP, schur.getP().getArray(), 1e-12);
    }

    @Test
    public void testVerifySchurDecomposition() {
        double[][] matrix = {
                {65, 35, 40, 69},
                {99, 64, 37, 2},
                {39, 48, 35, 90},
                {30, 93, 87, 17}
        };

        MatrixDense A = MatrixDense.from2DArray(matrix);
        SchurDecompositionDense schur = new SchurDecompositionDense(A);
        MatrixDense actual = schur.getP().multiply(schur.getT()).multiply(schur.getPT());

        assertArrayEquals(A.getArray(), actual.getArray(), 1e-12);
    }

    @Test
    public void allTests() {
        MatrixDense A, B, C, Z, O, I, R, S, X, SUB, M, T, SQ, DEF, SOL;
        // Uncomment this to test IO in a different locale.
        // Locale.setDefault(Locale.GERMAN);
        int errorCount = 0;
        int warningCount = 0;
        double tmp;
        double[] columnwise = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
        double[] rowwise = {1., 4., 7., 10., 2., 5., 8., 11., 3., 6., 9., 12.};
        double[][] avals = {{1., 4., 7., 10.}, {2., 5., 8., 11.}, {3., 6., 9., 12.}};
        double[][] rankdef = avals;
        double[][] tvals = {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}, {10., 11., 12.}};
        double[][] subavals = {{5., 8., 11.}, {6., 9., 12.}};
        double[][] pvals = {{4., 1., 1.}, {1., 2., 3.}, {1., 3., 6.}};
        double[][] ivals = {{1., 0., 0., 0.}, {0., 1., 0., 0.}, {0., 0., 1., 0.}};
        double[][] evals = {{0., 1., 0., 0.}, {1., 0., 2.e-7, 0.}, {0., -2.e-7, 0., 1.}, {0., 0., 1., 0.}};
        double[][] square = {{166., 188., 210.}, {188., 214., 240.}, {210., 240., 270.}};
        double[][] sqSolution = {{13.}, {15.}};
        double[][] condmat = {{1., 3.}, {7., 9.}};
        int rows = 3, cols = 4;
        int invalidld = 5;/* should trigger bad shape for construction with val */
        int validld = 3; /* leading dimension of intended test Matrices */
        int nonconformld = 4; /* leading dimension which is valid, but nonconforming */
        int ib = 1, ie = 2, jb = 1, je = 3; /* index ranges for sub Matrix */
        int[] rowindexset = {1, 2};
        int[] badrowindexset = {1, 3};
        int[] columnindexset = {1, 2, 3};
        int[] badcolumnindexset = {1, 2, 4};
        double columnsummax = 33.;
        double rowsummax = 30.;
        double sumofdiagonals = 15;
        double sumofsquares = 650;

        /**
         * Constructors and constructor-like methods: double[], int double[][] int, int
         * int, int, double int, int, double[][] constructWithCopy(double[][])
         * Random(int,int) identity(int)
         **/

        print("\nTesting constructors and constructor-like methods...\n");
        try {
            /** check that exception is thrown in packed constructor with invalid length **/
            new MatrixDense(columnwise, invalidld);
            errorCount = try_failure(errorCount, "Catch invalid length in packed constructor... ",
                    "exception not thrown for invalid input");
        } catch (IllegalArgumentException e) {
            try_success("Catch invalid length in packed constructor... ", e.getMessage());
        }
        try {
            /**
             * check that exception is thrown in default constructor if input array is
             * 'ragged'
             **/
        } catch (IllegalArgumentException e) {
            try_success("Catch ragged input to default constructor... ", e.getMessage());
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            errorCount = try_failure(errorCount, "Catch ragged input to constructor... ",
                    "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later");
        }
        try {
            /**
             * check that exception is thrown in constructWithCopy if input array is
             * 'ragged'
             **/
        } catch (IllegalArgumentException e) {
            try_success("Catch ragged input to constructWithCopy... ", e.getMessage());
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            errorCount = try_failure(errorCount, "Catch ragged input to constructWithCopy... ",
                    "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later");
        }

        A = new MatrixDense(columnwise, validld);
        B = new MatrixDense(avals);
        tmp = B.get(0, 0);
        avals[0][0] = 0.0;
        B.subtract(A);
        avals[0][0] = tmp;
        B = new MatrixDense(avals);
        tmp = B.get(0, 0);
        avals[0][0] = 0.0;
        if ((tmp - B.get(0, 0)) != 0.0) {
            /** check that constructWithCopy behaves properly **/
            errorCount = try_failure(errorCount, "constructWithCopy... ", "copy not effected... data visible outside");
        } else {
            try_success("constructWithCopy... ", "");
        }
        avals[0][0] = columnwise[0];
        I = new MatrixDense(ivals);
        try {
            check(I, MatrixDense.Factory.identity(3, 4));
            try_success("identity... ", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "identity... ", "identity Matrix not successfully created");
        }

        /**
         * Access Methods: getColumnCount() getRowCount() getArray() getArrayCopy()
         * getColumnPackedCopy() getRowPackedCopy() get(int,int)
         * subMatrix(int,int,int,int) subMatrix(int,int,int[]) subMatrix(int[],int,int)
         * subMatrix(int[],int[]) set(int,int,double) setMatrix(int,int,int,int,Matrix)
         * setMatrix(int,int,int[],Matrix) setMatrix(int[],int,int,Matrix)
         * setMatrix(int[],int[],Matrix)
         **/

        print("\nTesting access methods...\n");

        /**
         * Various get methods:
         **/

        B = new MatrixDense(avals);
        if (B.getRowCount() != rows) {
            errorCount = try_failure(errorCount, "getRowCount... ", "");
        } else {
            try_success("getRowCount... ", "");
        }
        if (B.getColumnCount() != cols) {
            errorCount = try_failure(errorCount, "getColumnCount... ", "");
        } else {
            try_success("getColumnCount... ", "");
        }
        B = new MatrixDense(avals);
        double[][] barray = B.getAs2DArray();
        if (barray == avals) {
            errorCount = try_failure(errorCount, "getArray... ", "");
        } else {
            try_success("getArray... ", "");
        }
        double[] bpacked = B.getColumnPackedCopy();
        try {
            check(bpacked, columnwise);
            try_success("getColumnPackedCopy... ", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "getColumnPackedCopy... ",
                    "data not successfully (deep) copied by columns");
        }
        bpacked = B.getRowPackedCopy();
        try {
            check(bpacked, rowwise);
            try_success("getRowPackedCopy... ", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "getRowPackedCopy... ", "data not successfully (deep) copied by rows");
        }
        try {
            B.get(B.getRowCount(), B.getColumnCount() - 1);
            errorCount = try_failure(errorCount, "get(int,int)... ", "OutOfBoundsException expected but not thrown");
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            try {
                B.get(B.getRowCount() - 1, B.getColumnCount());
                errorCount = try_failure(errorCount, "get(int,int)... ",
                        "OutOfBoundsException expected but not thrown");
            } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
                try_success("get(int,int)... OutofBoundsException... ", "");
            }
        } catch (java.lang.IllegalArgumentException e1) {
            errorCount = try_failure(errorCount, "get(int,int)... ", "OutOfBoundsException expected but not thrown");
        }
        try {
            if (B.get(B.getRowCount() - 1,
                    B.getColumnCount() - 1) != avals[B.getRowCount() - 1][B.getColumnCount() - 1]) {
                errorCount = try_failure(errorCount, "get(int,int)... ",
                        "Matrix entry (i,j) not successfully retrieved");
            } else {
                try_success("get(int,int)... ", "");
            }
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            errorCount = try_failure(errorCount, "get(int,int)... ", "Unexpected ArrayIndexOutOfBoundsException");
        }
        SUB = new MatrixDense(subavals);
        try {
            M = B.subMatrix(ib, ie + B.getRowCount() + 1, jb, je);
            errorCount = try_failure(errorCount, "subMatrix(int,int,int,int)... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            try {
                M = B.subMatrix(ib, ie, jb, je + B.getColumnCount() + 1);
                errorCount = try_failure(errorCount, "subMatrix(int,int,int,int)... ",
                        "ArrayIndexOutOfBoundsException expected but not thrown");
            } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
                try_success("subMatrix(int,int,int,int)... ArrayIndexOutOfBoundsException... ", "");
            }
        } catch (java.lang.IllegalArgumentException e1) {
            errorCount = try_failure(errorCount, "subMatrix(int,int,int,int)... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        }
        try {
            M = B.subMatrix(ib, ie, jb, je);
            try {
                check(SUB, M);
                try_success("subMatrix(int,int,int,int)... ", "");
            } catch (java.lang.RuntimeException e) {
                errorCount = try_failure(errorCount, "subMatrix(int,int,int,int)... ",
                        "sub-matrix not successfully retrieved");
            }
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            errorCount = try_failure(errorCount, "subMatrix(int,int,int,int)... ",
                    "Unexpected ArrayIndexOutOfBoundsException");
        }

        try {
            M = B.subMatrix(ib, ie, badcolumnindexset);
            errorCount = try_failure(errorCount, "subMatrix(int,int,int[])... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            try {
                M = B.subMatrix(ib, ie + B.getRowCount() + 1, columnindexset);
                errorCount = try_failure(errorCount, "subMatrix(int,int,int[])... ",
                        "ArrayIndexOutOfBoundsException expected but not thrown");
            } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
                try_success("subMatrix(int,int,int[])... ArrayIndexOutOfBoundsException... ", "");
            }
        } catch (java.lang.IllegalArgumentException e1) {
            errorCount = try_failure(errorCount, "subMatrix(int,int,int[])... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        }
        try {
            M = B.subMatrix(ib, ie, columnindexset);
            try {
                check(SUB, M);
                try_success("subMatrix(int,int,int[])... ", "");
            } catch (java.lang.RuntimeException e) {
                errorCount = try_failure(errorCount, "subMatrix(int,int,int[])... ",
                        "sub-matrix not successfully retrieved");
            }
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            errorCount = try_failure(errorCount, "subMatrix(int,int,int[])... ",
                    "Unexpected ArrayIndexOutOfBoundsException");
        }
        try {
            M = B.subMatrix(badrowindexset, jb, je);
            errorCount = try_failure(errorCount, "subMatrix(int[],int,int)... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            try {
                M = B.subMatrix(rowindexset, jb, je + B.getColumnCount() + 1);
                errorCount = try_failure(errorCount, "subMatrix(int[],int,int)... ",
                        "ArrayIndexOutOfBoundsException expected but not thrown");
            } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
                try_success("subMatrix(int[],int,int)... ArrayIndexOutOfBoundsException... ", "");
            }
        } catch (java.lang.IllegalArgumentException e1) {
            errorCount = try_failure(errorCount, "subMatrix(int[],int,int)... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        }
        try {
            M = B.subMatrix(rowindexset, jb, je);
            try {
                check(SUB, M);
                try_success("subMatrix(int[],int,int)... ", "");
            } catch (java.lang.RuntimeException e) {
                errorCount = try_failure(errorCount, "subMatrix(int[],int,int)... ",
                        "sub-matrix not successfully retrieved");
            }
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            errorCount = try_failure(errorCount, "subMatrix(int[],int,int)... ",
                    "Unexpected ArrayIndexOutOfBoundsException");
        }
        try {
            M = B.subMatrix(badrowindexset, columnindexset);
            errorCount = try_failure(errorCount, "subMatrix(int[],int[])... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            try {
                M = B.subMatrix(rowindexset, badcolumnindexset);
                errorCount = try_failure(errorCount, "subMatrix(int[],int[])... ",
                        "ArrayIndexOutOfBoundsException expected but not thrown");
            } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
                try_success("subMatrix(int[],int[])... ArrayIndexOutOfBoundsException... ", "");
            }
        } catch (java.lang.IllegalArgumentException e1) {
            errorCount = try_failure(errorCount, "subMatrix(int[],int[])... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        }
        try {
            M = B.subMatrix(rowindexset, columnindexset);
            try {
                check(SUB, M);
                try_success("subMatrix(int[],int[])... ", "");
            } catch (java.lang.RuntimeException e) {
                errorCount = try_failure(errorCount, "subMatrix(int[],int[])... ",
                        "sub-matrix not successfully retrieved");
            }
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            errorCount = try_failure(errorCount, "subMatrix(int[],int[])... ",
                    "Unexpected ArrayIndexOutOfBoundsException");
        }

        /**
         * Various set methods:
         **/

        try {
            B.set(B.getRowCount(), B.getColumnCount() - 1, 0.);
            errorCount = try_failure(errorCount, "set(int,int,double)... ",
                    "OutOfBoundsException expected but not thrown");
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            try {
                B.set(B.getRowCount() - 1, B.getColumnCount(), 0.);
                errorCount = try_failure(errorCount, "set(int,int,double)... ",
                        "OutOfBoundsException expected but not thrown");
            } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
                try_success("set(int,int,double)... OutofBoundsException... ", "");
            }
        } catch (java.lang.IllegalArgumentException e1) {
            errorCount = try_failure(errorCount, "set(int,int,double)... ",
                    "OutOfBoundsException expected but not thrown");
        }
        try {
            B.set(ib, jb, 0.);
            tmp = B.get(ib, jb);
            try {
                check(tmp, 0.);
                try_success("set(int,int,double)... ", "");
            } catch (java.lang.RuntimeException e) {
                errorCount = try_failure(errorCount, "set(int,int,double)... ", "Matrix element not successfully set");
            }
        } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
            errorCount = try_failure(errorCount, "set(int,int,double)... ",
                    "Unexpected ArrayIndexOutOfBoundsException");
        }
        M = new MatrixDense(2, 3, 0.);
        try {
            B.setMatrix(ib, ie + B.getRowCount() + 1, jb, je, M);
            errorCount = try_failure(errorCount, "setMatrix(int,int,int,int,Matrix)... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            try {
                B.setMatrix(ib, ie, jb, je + B.getColumnCount() + 1, M);
                errorCount = try_failure(errorCount, "setMatrix(int,int,int,int,Matrix)... ",
                        "ArrayIndexOutOfBoundsException expected but not thrown");
            } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
                try_success("setMatrix(int,int,int,int,Matrix)... ArrayIndexOutOfBoundsException... ", "");
            }
        } catch (java.lang.IllegalArgumentException e1) {
            errorCount = try_failure(errorCount, "setMatrix(int,int,int,int,Matrix)... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        }
        try {
            B.setMatrix(ib, ie, jb, je, M);
            try {
                check(M.subtract(B.subMatrix(ib, ie, jb, je)), M);
                try_success("setMatrix(int,int,int,int,Matrix)... ", "");
            } catch (java.lang.RuntimeException e) {
                errorCount = try_failure(errorCount, "setMatrix(int,int,int,int,Matrix)... ",
                        "sub-matrix not successfully set");
            }
            B.setMatrix(ib, ie, jb, je, SUB);
        } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
            errorCount = try_failure(errorCount, "setMatrix(int,int,int,int,Matrix)... ",
                    "Unexpected ArrayIndexOutOfBoundsException");
        }
        try {
            B.setMatrix(ib, ie + B.getRowCount() + 1, columnindexset, M);
            errorCount = try_failure(errorCount, "setMatrix(int,int,int[],Matrix)... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            try {
                B.setMatrix(ib, ie, badcolumnindexset, M);
                errorCount = try_failure(errorCount, "setMatrix(int,int,int[],Matrix)... ",
                        "ArrayIndexOutOfBoundsException expected but not thrown");
            } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
                try_success("setMatrix(int,int,int[],Matrix)... ArrayIndexOutOfBoundsException... ", "");
            }
        } catch (java.lang.IllegalArgumentException e1) {
            errorCount = try_failure(errorCount, "setMatrix(int,int,int[],Matrix)... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        }
        try {
            B.setMatrix(ib, ie, columnindexset, M);
            try {
                check(M.subtract(B.subMatrix(ib, ie, columnindexset)), M);
                try_success("setMatrix(int,int,int[],Matrix)... ", "");
            } catch (java.lang.RuntimeException e) {
                errorCount = try_failure(errorCount, "setMatrix(int,int,int[],Matrix)... ",
                        "sub-matrix not successfully set");
            }
            B.setMatrix(ib, ie, jb, je, SUB);
        } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
            errorCount = try_failure(errorCount, "setMatrix(int,int,int[],Matrix)... ",
                    "Unexpected ArrayIndexOutOfBoundsException");
        }
        try {
            B.setMatrix(rowindexset, jb, je + B.getColumnCount() + 1, M);
            errorCount = try_failure(errorCount, "setMatrix(int[],int,int,Matrix)... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            try {
                B.setMatrix(badrowindexset, jb, je, M);
                errorCount = try_failure(errorCount, "setMatrix(int[],int,int,Matrix)... ",
                        "ArrayIndexOutOfBoundsException expected but not thrown");
            } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
                try_success("setMatrix(int[],int,int,Matrix)... ArrayIndexOutOfBoundsException... ", "");
            }
        } catch (java.lang.IllegalArgumentException e1) {
            errorCount = try_failure(errorCount, "setMatrix(int[],int,int,Matrix)... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        }
        try {
            B.setMatrix(rowindexset, jb, je, M);
            try {
                check(M.subtract(B.subMatrix(rowindexset, jb, je)), M);
                try_success("setMatrix(int[],int,int,Matrix)... ", "");
            } catch (java.lang.RuntimeException e) {
                errorCount = try_failure(errorCount, "setMatrix(int[],int,int,Matrix)... ",
                        "sub-matrix not successfully set");
            }
            B.setMatrix(ib, ie, jb, je, SUB);
        } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
            errorCount = try_failure(errorCount, "setMatrix(int[],int,int,Matrix)... ",
                    "Unexpected ArrayIndexOutOfBoundsException");
        }
        try {
            B.setMatrix(rowindexset, badcolumnindexset, M);
            errorCount = try_failure(errorCount, "setMatrix(int[],int[],Matrix)... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            try {
                B.setMatrix(badrowindexset, columnindexset, M);
                errorCount = try_failure(errorCount, "setMatrix(int[],int[],Matrix)... ",
                        "ArrayIndexOutOfBoundsException expected but not thrown");
            } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
                try_success("setMatrix(int[],int[],Matrix)... ArrayIndexOutOfBoundsException... ", "");
            }
        } catch (java.lang.IllegalArgumentException e1) {
            errorCount = try_failure(errorCount, "setMatrix(int[],int[],Matrix)... ",
                    "ArrayIndexOutOfBoundsException expected but not thrown");
        }
        try {
            B.setMatrix(rowindexset, columnindexset, M);
            try {
                check(M.subtract(B.subMatrix(rowindexset, columnindexset)), M);
                try_success("setMatrix(int[],int[],Matrix)... ", "");
            } catch (java.lang.RuntimeException e) {
                errorCount = try_failure(errorCount, "setMatrix(int[],int[],Matrix)... ",
                        "sub-matrix not successfully set");
            }
        } catch (java.lang.ArrayIndexOutOfBoundsException e1) {
            errorCount = try_failure(errorCount, "setMatrix(int[],int[],Matrix)... ",
                    "Unexpected ArrayIndexOutOfBoundsException");
        }

        /**
         * Array-like methods: subtract subtractEquals plus addEquals arrayLeftDivide
         * arrayLeftDivideEquals arrayRightDivide arrayRightDivideEquals arrayMultiply
         * arrayMultiplyEquals usubtract
         **/

        print("\nTesting array-like methods...\n");
        S = new MatrixDense(columnwise, nonconformld);
        R = MatrixDense.Factory.random(A.getRowCount(), A.getColumnCount());
        A = R;
        try {
            S = A.subtract(S);
            errorCount = try_failure(errorCount, "subtract conformance check... ", "nonconformance not raised");
        } catch (IllegalArgumentException e) {
            try_success("subtract conformance check... ", "");
        }
        if (A.subtract(R).norm1() != 0.) {
            errorCount = try_failure(errorCount, "subtract... ",
                    "(difference of identical Matrices is nonzero,\nSubsequent use of subtract should be suspect)");
        } else {
            try_success("subtract... ", "");
        }
        A = R.copy();
        A.subtractEquals(R);
        Z = new MatrixDense(A.getRowCount(), A.getColumnCount());
        try {
            A.subtractEquals(S);
            errorCount = try_failure(errorCount, "subtractEquals conformance check... ", "nonconformance not raised");
        } catch (IllegalArgumentException e) {
            try_success("subtractEquals conformance check... ", "");
        }
        if (A.subtract(Z).norm1() != 0.) {
            errorCount = try_failure(errorCount, "subtractEquals... ",
                    "(difference of identical Matrices is nonzero,\nSubsequent use of subtract should be suspect)");
        } else {
            try_success("subtractEquals... ", "");
        }

        A = R.copy();
        B = MatrixDense.Factory.random(A.getRowCount(), A.getColumnCount());
        C = A.subtract(B);
        try {
            S = A.add(S);
            errorCount = try_failure(errorCount, "plus conformance check... ", "nonconformance not raised");
        } catch (IllegalArgumentException e) {
            try_success("plus conformance check... ", "");
        }
        try {
            check(C.add(B), A);
            try_success("plus... ", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "plus... ", "(C = A - B, but C + B != A)");
        }
        C = A.subtract(B);
        C.addEquals(B);
        try {
            A.addEquals(S);
            errorCount = try_failure(errorCount, "addEquals conformance check... ", "nonconformance not raised");
        } catch (IllegalArgumentException e) {
            try_success("addEquals conformance check... ", "");
        }
        try {
            check(C, A);
            try_success("addEquals... ", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "addEquals... ", "(C = A - B, but C = C + B != A)");
        }
        A = R.uminus();
        try {
            check(A.add(R), Z);
            try_success("usubtract... ", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "usubtract... ", "(-A + A != zeros)");
        }
        A = R.copy();
        O = new MatrixDense(A.getRowCount(), A.getColumnCount(), 1.0);
        C = A.arrayLeftDivide(R);
        try {
            S = A.arrayLeftDivide(S);
            errorCount = try_failure(errorCount, "arrayLeftDivide conformance check... ", "nonconformance not raised");
        } catch (IllegalArgumentException e) {
            try_success("arrayLeftDivide conformance check... ", "");
        }
        try {
            check(C, O);
            try_success("arrayLeftDivide... ", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "arrayLeftDivide... ", "(M.\\M != ones)");
        }
        try {
            A.arrayLeftDivideEquals(S);
            errorCount = try_failure(errorCount, "arrayLeftDivideEquals conformance check... ",
                    "nonconformance not raised");
        } catch (IllegalArgumentException e) {
            try_success("arrayLeftDivideEquals conformance check... ", "");
        }
        A.arrayLeftDivideEquals(R);
        try {
            check(A, O);
            try_success("arrayLeftDivideEquals... ", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "arrayLeftDivideEquals... ", "(M.\\M != ones)");
        }
        A = R.copy();
        try {
            A.arrayRightDivide(S);
            errorCount = try_failure(errorCount, "arrayRightDivide conformance check... ", "nonconformance not raised");
        } catch (IllegalArgumentException e) {
            try_success("arrayRightDivide conformance check... ", "");
        }
        C = A.arrayRightDivide(R);
        try {
            check(C, O);
            try_success("arrayRightDivide... ", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "arrayRightDivide... ", "(M./M != ones)");
        }
        try {
            A.arrayRightDivideEquals(S);
            errorCount = try_failure(errorCount, "arrayRightDivideEquals conformance check... ",
                    "nonconformance not raised");
        } catch (IllegalArgumentException e) {
            try_success("arrayRightDivideEquals conformance check... ", "");
        }
        A.arrayRightDivideEquals(R);
        try {
            check(A, O);
            try_success("arrayRightDivideEquals... ", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "arrayRightDivideEquals... ", "(M./M != ones)");
        }
        A = R.copy();
        B = MatrixDense.Factory.random(A.getRowCount(), A.getColumnCount());
        try {
            S = A.arrayMultiply(S);
            errorCount = try_failure(errorCount, "arrayMultiply conformance check... ", "nonconformance not raised");
        } catch (IllegalArgumentException e) {
            try_success("arrayMultiply conformance check... ", "");
        }
        C = A.arrayMultiply(B);
        try {
            C.arrayRightDivideEquals(B);
            check(C, A);
            try_success("arrayMultiply... ", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "arrayMultiply... ", "(A = R, C = A.*B, but C./B != A)");
        }
        try {
            A.arrayMultiplyEquals(S);
            errorCount = try_failure(errorCount, "arrayMultiplyEquals conformance check... ",
                    "nonconformance not raised");
        } catch (IllegalArgumentException e) {
            try_success("arrayMultiplyEquals conformance check... ", "");
        }
        A.arrayMultiplyEquals(B);
        try {
            A.arrayRightDivideEquals(B);
            check(A, R);
            try_success("arrayMultiplyEquals... ", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "arrayMultiplyEquals... ", "(A = R, A = A.*B, but A./B != R)");
        }

        /**
         * LA methods: transpose multiply cond rank det trace norm1 norm2 normF normInf
         * solve solveTranspose inverse chol eig lu qr svd
         **/

        print("\nTesting linear algebra methods...\n");
        A = new MatrixDense(columnwise, 3);
        T = new MatrixDense(tvals);
        T = A.transpose();
        try {
            check(A.transpose(), T);
            try_success("transpose...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "transpose()...", "transpose unsuccessful");
        }
        A.transpose();
        try {
            check(A.norm1(), columnsummax);
            try_success("norm1...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "norm1()...", "incorrect norm calculation");
        }
        try {
            check(A.normInf(), rowsummax);
            try_success("normInf()...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "normInf()...", "incorrect norm calculation");
        }
        try {
            check(A.normFrob(), Math.sqrt(sumofsquares));
            try_success("normF...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "normF()...", "incorrect norm calculation");
        }
        try {
            check(A.trace(), sumofdiagonals);
            try_success("trace()...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "trace()...", "incorrect trace calculation");
        }
        try {
            check(A.subMatrix(0, A.getRowCount() - 1, 0, A.getRowCount() - 1).det(), 0.);
            try_success("det()...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "det()...", "incorrect determinant calculation");
        }
        SQ = new MatrixDense(square);
        try {
            check(A.multiply(A.transpose()), SQ);
            try_success("multiply(Matrix)...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "multiply(Matrix)...", "incorrect Matrix-Matrix product calculation");
        }
        try {
            check(A.multiply(0.), Z);
            try_success("multiply(double)...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "multiply(double)...", "incorrect Matrix-scalar product calculation");
        }

        A = new MatrixDense(columnwise, 4);
        QRDecompositionDense QR = A.QR();
        R = QR.getR();
        try {
            check(A, QR.getQThin().multiply(R));
            try_success("QRDecomposition...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "QRDecomposition...", "incorrect QR decomposition calculation");
        }
        SingularValueDecompositionDense SVD = A.SVD();
        try {
            check(A, SVD.getU().multiply(SVD.getS().multiply(SVD.getV().transpose())));
            try_success("SingularValueDecomposition...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "SingularValueDecomposition...",
                    "incorrect singular value decomposition calculation");
        }
        DEF = new MatrixDense(rankdef);
        try {
            check(DEF.rank(), Math.min(DEF.getRowCount(), DEF.getColumnCount()) - 1);
            try_success("rank()...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "rank()...", "incorrect rank calculation");
        }
        B = new MatrixDense(condmat);
        SVD = B.SVD();
        double[] singularvalues = SVD.getSingularValues();
        try {
            check(B.cond(), singularvalues[0] / singularvalues[Math.min(B.getRowCount(), B.getColumnCount()) - 1]);
            try_success("cond()...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "cond()...", "incorrect condition number calculation");
        }
        int n = A.getColumnCount();
        A = A.subMatrix(0, n - 1, 0, n - 1);
        A.set(0, 0, 0.);
        LUDecompositionDense LU = A.LU();
        try {
            check(A.subMatrix(LU.getPivot(), 0, n - 1), LU.getL().multiply(LU.getU()));
            try_success("LUDecomposition...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "LUDecomposition...", "incorrect LU decomposition calculation");
        }
        X = A.inv();
        try {
            check(A.multiply(X), MatrixDense.Factory.identity(3));
            try_success("inverse()...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "inverse()...", "incorrect inverse calculation");
        }
        O = new MatrixDense(SUB.getRowCount(), 1, 1.0);
        SOL = new MatrixDense(sqSolution);
        SQ = SUB.subMatrix(0, SUB.getRowCount() - 1, 0, SUB.getRowCount() - 1);
        try {
            check(SQ.solve(SOL), O);
            try_success("solve()...", "");
        } catch (java.lang.IllegalArgumentException e1) {
            errorCount = try_failure(errorCount, "solve()...", e1.getMessage());
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "solve()...", e.getMessage());
        }
        A = new MatrixDense(pvals);
        CholeskyDecompositionDense Chol = A.Chol();
        MatrixDense L = Chol.getL();
        try {
            check(A, L.multiply(L.transpose()));
            try_success("CholeskyDecomposition...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "CholeskyDecomposition...",
                    "incorrect Cholesky decomposition calculation");
        }
        X = Chol.solve(MatrixDense.Factory.identity(3));
        try {
            check(A.multiply(X), MatrixDense.Factory.identity(3));
            try_success("CholeskyDecomposition solve()...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "CholeskyDecomposition solve()...",
                    "incorrect Choleskydecomposition solve calculation");
        }
        EigenvalueDecompositionDense Eig = A.eig();
        MatrixDense D = Eig.getD();
        MatrixDense V = Eig.getV();
        try {
            check(A.multiply(V), V.multiply(D));
            try_success("EigenvalueDecomposition (symmetric)...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "EigenvalueDecomposition (symmetric)...",
                    "incorrect symmetric Eigenvalue decomposition calculation");
        }
        A = new MatrixDense(evals);
        Eig = A.eig();
        D = Eig.getD();
        V = Eig.getV();
        try {
            check(A.multiply(V), V.multiply(D));
            try_success("EigenvalueDecomposition (nonsymmetric)...", "");
        } catch (java.lang.RuntimeException e) {
            errorCount = try_failure(errorCount, "EigenvalueDecomposition (nonsymmetric)...",
                    "incorrect nonsymmetric Eigenvalue decomposition calculation");
        }

        assertEquals(0, errorCount);

        print("\nTestMatrix completed.\n");
        print("Total errors reported: " + Integer.toString(errorCount) + "\n");
        print("Total warnings reported: " + Integer.toString(warningCount) + "\n");
    }

    /** private utility routines **/

    /**
     * Check magnitude of difference of scalars.
     **/

    private static void check(double x, double y) {
        double eps = ConstantsETK.DOUBLE_EPS;
        if (x == 0 & Math.abs(y) < 10 * eps)
            return;
        if (y == 0 & Math.abs(x) < 10 * eps)
            return;
        if (Math.abs(x - y) > 10 * eps * Math.max(Math.abs(x), Math.abs(y))) {
            throw new RuntimeException(
                    "The difference x-y is too large: x = " + Double.toString(x) + "  y = " + Double.toString(y));
        }
    }

    /**
     * Check norm of difference of "vectors".
     **/

    private static void check(double[] x, double[] y) {
        if (x.length == y.length) {
            for (int i = 0; i < x.length; i++) {
                check(x[i], y[i]);
            }
        } else {
            throw new RuntimeException("Attempt to compare vectors of different lengths");
        }
    }

    /**
     * Check norm of difference of Matrices.
     **/

    private static void check(MatrixDense X, MatrixDense Y) {
        double eps = ConstantsETK.DOUBLE_EPS;
        if (X.norm1() == 0. & Y.norm1() < 10 * eps)
            return;
        if (Y.norm1() == 0. & X.norm1() < 10 * eps)
            return;
        if (X.subtract(Y).norm1() > 1000 * eps * Math.max(X.norm1(), Y.norm1())) {
            throw new RuntimeException("The norm of (X-Y) is too large: " + Double.toString(X.subtract(Y).norm1()));
        }
    }

    /**
     * Shorten spelling of print.
     **/

    private static void print(String s) {
//        Enable to see test output
//        System.out.print(s);
    }

    /**
     * Print appropriate messages for successful outcome try
     **/

    private static void try_success(String s, String e) {
        print(">    " + s + "success\n");
        if (e != "") {
            print(">      Message: " + e + "\n");
        }
    }

    /**
     * Print appropriate messages for unsuccessful outcome try
     **/

    private static int try_failure(int count, String s, String e) {
        print(">    " + s + "*** failure ***\n>      Message: " + e + "\n");
        return ++count;
    }
}
