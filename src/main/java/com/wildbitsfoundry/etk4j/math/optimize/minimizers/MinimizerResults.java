package com.wildbitsfoundry.etk4j.math.optimize.minimizers;

import com.wildbitsfoundry.etk4j.math.optimize.OptimizerStatusType;

public class MinimizerResults<T> {
    private T value;
    private T functionValue;
    private int numberOfIterations;
    private String minimizerStatus;
    private boolean converged;
    private OptimizerStatusType optimizerStatusType;

    public T getValue() {
        return value;
    }

    void setValue(T value) {
        this.value = value;
    }

    public T getFunctionValue() {
        return functionValue;
    }

    public void setFunctionValue(T functionValue) {
        this.functionValue = functionValue;
    }

    public int getNumberOfIterations() {
        return numberOfIterations;
    }

    void setNumberOfIterations(int numberOfIterations) {
        this.numberOfIterations = numberOfIterations;
    }

    public String getMinimizerStatus() {
        return minimizerStatus;
    }

    void setMinimizerStatus(String minimizerStatus) {
        this.minimizerStatus = minimizerStatus;
    }

    public OptimizerStatusType getOptimizerStatusType() {
        return optimizerStatusType;
    }

    public void setOptimizerStatusType(OptimizerStatusType optimizerStatusType) {
        this.optimizerStatusType = optimizerStatusType;
    }

    public boolean hasConverged() {
        return converged;
    }

    void setHasConverged(boolean converged) {
        this.converged = converged;
    }

    @Override
    public String toString() {
        return "SolverResults{" +
                "value=" + value +
                ", functionValue=" + functionValue +
                ", numberOfIterations=" + numberOfIterations +
                ", solverStatus='" + minimizerStatus + '\'' +
                ", converged=" + converged +
                '}';
    }
}

