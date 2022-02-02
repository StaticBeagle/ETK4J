package com.wildbitsfoundry.etk4j.control;

public class NonUniformTimeStepsException extends RuntimeException {
    public NonUniformTimeStepsException(){}

    public NonUniformTimeStepsException(String message) {
        super(message);
    }

    public NonUniformTimeStepsException(String message, Throwable cause) {
        super(message, cause);
    }

    public NonUniformTimeStepsException(Throwable cause) {
        super(cause);
    }
}
