package com.wildbitsfoundry.etk4j.util.data;

public class DoubleDynamicArray {

    private int count;
    private int capacity;
    private double[] data;

    public DoubleDynamicArray() {
        this(10);
    }

    public DoubleDynamicArray(int capacity) {
        if(capacity < 0) {
            throw new IllegalArgumentException("Initial capacity must be greater than 0");
        }
        this.count = 0;
        this.capacity = capacity;
        this.data = new double[capacity];
    }

    public DoubleDynamicArray(double ...data) {
        this.count = data.length;
        this.capacity = data.length;
        System.arraycopy(data, 0, this.data, 0, data.length);
    }

    public int size() {
        return count;
    }

    public void add(double d) {
        if (count++ == capacity) {
            if(capacity == 0) {
                capacity++;
            }
            capacity += capacity >> 2;
            growContainer(capacity);
        }
        data[count - 1] = d;
    }

    public double get(int index) {
        return data[index];
    }

    public void set(int index, double val) {
        if (index >= capacity) {
            throw new IllegalArgumentException("Out of bounds"); // better exception add
        }
        data[index] = val;
    }

    public void clear() {
        this.count = 0;
        this.capacity = 10;
        this.data = new double[10];
    }

    public double[] getArray() {
        return data;
    }

    public void reshape(int length) {
        if (data.length < length) {
            data = new double[length];
        }
        this.count = length;
    }

    private void growContainer(int size) {
        double[] data = new double[size];
        System.arraycopy(this.data, 0, data, 0, this.data.length);
        this.data = data;
    }
}
