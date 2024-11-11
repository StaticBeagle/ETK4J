package com.wildbitsfoundry.etk4j.util.data;

public class IntDynamicArray {

    private int count;
    private int capacity;
    private int[] data;

    public IntDynamicArray() {
        this(10);
    }

    public IntDynamicArray(int capacity) {
        this.count = 0;
        this.capacity = capacity;
        this.data = new int[capacity];
    }

    public int size() {
        return count;
    }

    public void add(int i) {
        if(count++ == capacity) {
            capacity += capacity >> 2;
            growContainer(capacity);
        }
        data[count - 1] = i;
    }

    public double get(int index) {
        return data[index];
    }

    public void set(int index, int val) {
        if(index >= capacity) {
            throw new IllegalArgumentException("Out of bounds"); // better exception add
        }
        data[index] = val;
    }

    public int[] getArray() {
        return data;
    }

    public void reshape(int length) {
        if (data.length < length) {
            data = new int[length];
        }
        this.count = length;
    }

    private void growContainer(int size) {
        int[] data = new int[size];
        System.arraycopy(this.data, 0, data, 0, this.data.length);
        this.data = data;
    }
}
