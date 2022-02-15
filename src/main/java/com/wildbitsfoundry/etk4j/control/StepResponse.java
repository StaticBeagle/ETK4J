package com.wildbitsfoundry.etk4j.control;

/**
 * The {@code StepResponse} class holds the time domain results of a continuous-time system.
 */
public class StepResponse {

    private double[] time;
    private double[] response;

    StepResponse(double[] time, double[] response) {
        this.time = time;
        this.response = response;
    }

    /**
     * Array of time points.
     * @return The times at which the response was calculated.
     */
    public double[] getTime() {
        return time;
    }

    /**
     * Time domain response.
     * @return The time domain response of the system.
     */
    public double[] getResponse() {
        return response;
    }
}
