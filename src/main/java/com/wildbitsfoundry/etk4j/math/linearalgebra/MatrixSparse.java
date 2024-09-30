package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.util.data.DynamicArray;

import java.util.Arrays;

public class MatrixSparse {

    private DynamicArray nonZeroValues;
    private int[] nonZeroRows;
    private int[] columnIndex;

    private int rows;
    private int cols;

    private boolean indicesSorted = false;

    public MatrixSparse(int rows, int cols) {
        this(rows, cols, 10);
    }

    public MatrixSparse(int rows, int cols, int nonZeroValuesLength) {
        // validate inputs
        this.nonZeroValues = new DynamicArray(nonZeroValuesLength);
        this.nonZeroRows = new int[nonZeroValuesLength];
        this.rows = rows;
        this.cols = cols;
        this.columnIndex = new int[cols + 1];
    }

    public boolean isEmpty() {
        return this.nonZeroValues.size() == 0;
    }

    public int getRowCount() {
        return rows;
    }

    public int getColumnCount() {
        return cols;
    }

    public double get(int row, int col) {
        return get(row, col, 0.0);
    }

    public double get(int row, int col, double valueIfNotFound) {
        // validate
        int index = findIndexInNonZeroRows(row, col);
        return index >= 0 ? nonZeroValues.get(index) : valueIfNotFound;
    }

    public void set(int row, int col, double val) {
        int index = findIndexInNonZeroRows(row, col);
        if(index >= 0) {
            nonZeroValues.set(index, val);
        } else {
            int idx0 = columnIndex[col];
            int idx1 = columnIndex[col + 1];

            index = idx0;
            while(index < idx1 && row >= nonZeroRows[index]) {
                index++;
            }

            for(int i = col + 1; i <= this.cols; i++) {
                this.columnIndex[i]++;
            }

            nonZeroValues.add(0);
            for(int i = nonZeroValues.size(); i > index; i--) {
                nonZeroRows[i] = nonZeroRows[i - 1];
                nonZeroValues.set(i, nonZeroValues.get(i - 1));
            }

            nonZeroRows[index] = row;
            nonZeroValues.set(index, val);

        }

        // implement sortIndices
    }

    private int findIndexInNonZeroRows(int row, int col) {
        int col0 = columnIndex[col];
        int col1 = columnIndex[col + 1];

        if(indicesSorted) {
            return Arrays.binarySearch(nonZeroRows, col0, col1, row);
        }

        for(int i = col0; i < col1; i++) {
            if(nonZeroRows[i] == row) {
                return i;
            }
        }
        return -1;
    }

    @Override
    public String toString() {
        if(isEmpty()) {
            return "[]";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                int index = findIndexInNonZeroRows(i, j);
                if(index >= 0) {
                    sb.append(String.format("%.4g", nonZeroValues.get(index))).append(" ");
                } else {
                    sb.append(". ");
                }
            }
            sb.append(System.lineSeparator());
        }
        sb.setLength(sb.length() - 1);
        return sb.toString();
    }

    public static void main(String[] args) {
        MatrixSparse matrix = new MatrixSparse(3, 3);

        matrix.set(0, 0, 1.0);
        matrix.set(1, 1, 2.0);
        matrix.set(2, 2, 3.0);

        System.out.println(matrix);
    }
}
