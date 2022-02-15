package com.wildbitsfoundry.etk4j.math.calculus;

@SuppressWarnings("serial")
/**
 * The {@code MaximumNumberOfEvaluationsReached} class is an {@link RuntimeException} that represents the case where
 * the number of evaluations of an integration algorithm exceeds the maximum number of allowed evaluations.
 */
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
