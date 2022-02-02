package com.wildbitsfoundry.etk4j.math.calculus;

@SuppressWarnings("serial")
public class MaximumNumberOfEvaluationsReached extends RuntimeException {
    public MaximumNumberOfEvaluationsReached() {
    }

    public MaximumNumberOfEvaluationsReached(String message) {
        super(message);
    }

    public MaximumNumberOfEvaluationsReached(String message, Throwable cause) {
        super(message, cause);
    }

    public MaximumNumberOfEvaluationsReached(Throwable cause) {
        super(cause);
    }

}
