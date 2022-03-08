package com.wildbitsfoundry.etk4j.signals.filters;

/**
 * The {@code NonFiniteRootsException} class is thrown when a collection of complex values contain values that are not
 * finite.
 */
public class NonFiniteRootsException extends RuntimeException {
    public NonFiniteRootsException(){}

    public NonFiniteRootsException(String message) {
        super(message);
    }

    public NonFiniteRootsException(String message, Throwable cause) {
        super(message, cause);
    }

    public NonFiniteRootsException(Throwable cause) {
        super(cause);
    }
}
