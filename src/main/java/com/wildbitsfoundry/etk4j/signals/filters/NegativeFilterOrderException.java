package com.wildbitsfoundry.etk4j.signals.filters;

// TODO other methods to match other exceptions
public class NegativeFilterOrderException extends RuntimeException {

    /**
	 * 
	 */
	private static final long serialVersionUID = 3135934681615853289L;

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
