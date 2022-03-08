package com.wildbitsfoundry.etk4j.signals.filters;

/**
 * The {@code NegativeFilterOrderException} class is thrown when a collection of arrays contain complex values that are not
 * finite.
 */
public class NegativeFilterOrderException extends RuntimeException {

    public NegativeFilterOrderException(){}

    public NegativeFilterOrderException(String message) {
        super(message);
    }

    public NegativeFilterOrderException(String message, Throwable cause) {
        super(message, cause);
    }

    public NegativeFilterOrderException(Throwable cause) {
        super(cause);
    }
}
