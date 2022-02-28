package com.wildbitsfoundry.etk4j.control;

/**
 * The {@code TimeResponse} class holds the time response of a continuous-time system.
 */
public class TimeResponse {

    private double[] time;
    private double[][] response;
    private double[][] evolutionOfStateVector;

    public TimeResponse(double[] time, double[][] response, double[][] evolutionOfStateVector) {
        this.time = time;
        this.response = response;
        this.evolutionOfStateVector = evolutionOfStateVector;
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
     * @return The time domain response of the system. Each row represents the time domain response for each output
     * of the multi-output system.
     */
    public double[][] getResponse() {
        return response;
    }

    /**
     * Evolution of the state vector
     * @return The evolution of the state vector vs time.
     */
    public double[][] getEvolutionOfStateVector() {
        return evolutionOfStateVector;
    }
}
