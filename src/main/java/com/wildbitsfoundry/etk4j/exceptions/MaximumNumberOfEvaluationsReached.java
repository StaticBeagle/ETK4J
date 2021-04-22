package com.wildbitsfoundry.etk4j.exceptions;

@SuppressWarnings("serial")
public class MaximumNumberOfEvaluationsReached extends RuntimeException {

    public MaximumNumberOfEvaluationsReached(String message) {
        super(message);
    }
}
