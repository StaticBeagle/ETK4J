package com.wildbitsfoundry.etk4j.control;

/**
 * The {@code ImproperTransferFunctionException} class is thrown when the number of the denominator is less than
 * the number of the numerator in a {@link TransferFunction} system. For many operations in a transfer function, it is
 * required that for the transfer function to be proper.
 */
public class ImproperTransferFunctionException extends RuntimeException {
    public ImproperTransferFunctionException(){}

    public ImproperTransferFunctionException(String message) {
        super(message);
    }

    public ImproperTransferFunctionException(String message, Throwable cause) {
        super(message, cause);
    }

    public ImproperTransferFunctionException(Throwable cause) {
        super(cause);
    }
}
