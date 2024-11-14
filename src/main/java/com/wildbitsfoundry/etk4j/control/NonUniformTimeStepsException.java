package com.wildbitsfoundry.etk4j.control;

/**
 * The {@code NonUniformTimeStepsException} is thrown when the time vector to a time simulation for a
 * {@link LinearTimeInvariantSystem} is not uniform i.e. the difference between each point is different.
 */
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
