package com.wildbitsfoundry.etk4j.signals.filters;

/**
 * The {@code MaximumNumberOfIterationsReachedException} class is thrown when the maximum number of iterations for an
 * iterative method is reached.
 */
public class MaximumNumberOfIterationsReachedException extends RuntimeException {
    public MaximumNumberOfIterationsReachedException(){}

    public MaximumNumberOfIterationsReachedException(String message) {
        super(message);
    }

    public MaximumNumberOfIterationsReachedException(String message, Throwable cause) {
        super(message, cause);
    }

    public MaximumNumberOfIterationsReachedException(Throwable cause) {
        super(cause);
    }
}
