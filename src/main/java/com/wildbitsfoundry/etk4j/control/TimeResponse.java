package com.wildbitsfoundry.etk4j.control;

public class TimeResponse {

    private double[] time;
    private double[][] response;
    private double[][] evolutionOfStateVector;

    public TimeResponse(double[] time, double[][] response, double[][] evolutionOfStateVector) {
        this.time = time;
        this.response = response;
        this.evolutionOfStateVector = evolutionOfStateVector;
    }

    public double[] getTime() {
        return time;
    }

    public double[][] getResponse() {
        return response;
    }

    public double[][] getEvolutionOfStateVector() {
        return evolutionOfStateVector;
    }
}
