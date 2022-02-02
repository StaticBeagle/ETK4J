package com.wildbitsfoundry.etk4j.control;

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
