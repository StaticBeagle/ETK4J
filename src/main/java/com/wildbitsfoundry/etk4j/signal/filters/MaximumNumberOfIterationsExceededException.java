package com.wildbitsfoundry.etk4j.signal.filters;

/**
 * The {@code MaximumNumberOfIterationsReachedException} class is thrown when the maximum number of iterations for an
 * iterative method is reached.
 */
public class MaximumNumberOfIterationsExceededException extends RuntimeException {
    public MaximumNumberOfIterationsExceededException(){}

    public MaximumNumberOfIterationsExceededException(String message) {
        super(message);
    }

    public MaximumNumberOfIterationsExceededException(String message, Throwable cause) {
        super(message, cause);
    }

    public MaximumNumberOfIterationsExceededException(Throwable cause) {
        super(cause);
    }
}
