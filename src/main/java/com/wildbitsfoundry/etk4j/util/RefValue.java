package com.wildbitsfoundry.etk4j.util;

public class RefValue<T> {
    private T value;

    public RefValue(T value) {
        this.value = value;
    }

    public void setValue(T value) {
        this.value = value;
    }

    public T getValue() {
        return this.value;
    }

    public static class RefDouble extends RefValue<Double> {
        public RefDouble(double value) {
            super(value);
        }
    }
    
    public static class RefInteger extends RefValue<Integer> {
    	public RefInteger(int value) {
    		super(value);
    	}
    }
}
