package com.wildbitsfoundry.etk4j.iterationtools;

import java.util.Arrays;
import java.util.Iterator;

public class CartesianProduct implements Iterable<int[]>, Iterator<int[]> {

    private final int[] _lengths;
    private final int[] _indices;
    private boolean _hasNext = true;

    public CartesianProduct(int[] lengths) {
        _lengths = lengths;
        _indices = new int[lengths.length];
    }

    public boolean hasNext() {
        return _hasNext;
    }

    public int[] next() {
        int[] result = Arrays.copyOf(_indices, _indices.length);
        for (int i = _indices.length - 1; i >= 0; i--) {
            if (_indices[i] == _lengths[i] - 1) {
                _indices[i] = 0;
                if (i == 0) {
                    _hasNext = false;
                }
            } else {
                _indices[i]++;
                break;
            }
        }
        return result;
    }

    public Iterator<int[]> iterator() {
        return this;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }
}