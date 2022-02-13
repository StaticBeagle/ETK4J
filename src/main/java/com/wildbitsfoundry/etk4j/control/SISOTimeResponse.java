package com.wildbitsfoundry.etk4j.control;

/**
 * The {@code SISOTimeResponse}  represents the time response af a Single Input Single Output system e.g.
 * {@link TransferFunction}.
 */
public class SISOTimeResponse {

    private double[] time;
    private double[] response;
    private double[][] evolutionOfStateVector;

    SISOTimeResponse(double[] time, double[] response, double[][] evolutionOfStateVector) {
        this.time = time;
        this.response = response;
        this.evolutionOfStateVector = evolutionOfStateVector;
    }

    /**
     * Time array.
     * @return The array of time points at which the response was evaluated.
     */
    public double[] getTime() {
        return time;
    }

    /**
     * Response vector.
     * @return The time response of the system.
     */
    public double[] getResponse() {
        return response;
    }

    /**
     * Evolution of State Vector.
     * @return The evolutions of all the values that the state vector went through while calculating the time response.
     */
    public double[][] getEvolutionOfStateVector() {
        return evolutionOfStateVector;
    }
}
