package com.wildbitsfoundry.etk4j.util.data;

public class DynamicArray {

    private int count;
    private int capacity;
    private double[] data;

    public DynamicArray() {
        this(10);
    }

    public DynamicArray(int capacity) {
        this.count = 0;
        this.capacity = capacity;
        this.data = new double[capacity];
    }

    public int size() {
        return count;
    }

    public void add(double d) {
        if(count++ == capacity) {
            capacity += capacity >> 2;
            growContainer(capacity);
        }
        data[count - 1] = d;
    }

    public double get(int index) {
        return data[index];
    }

    public void set(int index, double val) {
        if(index >= capacity) {
            throw new IllegalArgumentException("Out of bounds"); // better exception add
        }
        data[index] = val;
    }

    public double[] getArray() {
        return data;
    }

    private void growContainer(int size) {
        double[] data = new double[size];
        System.arraycopy(this.data, 0, data, 0, this.data.length);
        this.data = data;
    }
}
