package com.wildbitsfoundry.etk4j.control;

public class TimeResponseResults {

    private double[] time;
    private double[] response;
    private double[][] stateVector;

    public TimeResponseResults(double[] time, double[] response, double[][] stateVector) {
        this.time = time;
        this.response = response;
        this.stateVector = stateVector;
    }

    public double[] getTime() {
        return time;
    }

    public double[] getResponse() {
        return response;
    }

    public double[][] getStateVector() {
        return stateVector;
    }
}
