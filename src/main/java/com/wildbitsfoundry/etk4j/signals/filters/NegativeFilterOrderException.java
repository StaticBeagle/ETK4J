package com.wildbitsfoundry.etk4j.signals.filters;

public class NegativeFilterOrderException extends RuntimeException {

    public NegativeFilterOrderException(String message) {
        super(message);
    }

    public NegativeFilterOrderException(Throwable cause) {
        super(cause);
    }

    public NegativeFilterOrderException(String message, Throwable cause) {
        super(message, cause);
    }
}
