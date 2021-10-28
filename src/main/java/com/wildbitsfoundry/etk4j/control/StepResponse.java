package com.wildbitsfoundry.etk4j.control;

public class StepResponse {

    private double[] time;
    private double[] response;

    public StepResponse(double[] time, double[] response) {
        this.time = time;
        this.response = response;
    }

    public double[] getTime() {
        return time;
    }

    public double[] getResponse() {
        return response;
    }
}
