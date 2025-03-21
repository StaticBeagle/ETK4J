package com.wildbitsfoundry.etk4j.signal.smoothing;

/**
 * An &alpha;-&beta; filter (also called alpha-beta filter, f-g filter or g-h filter) is a simplified form of observer
 * for estimation, data smoothing and control applications. It is closely related to Kalman filters and to linear state
 * observers used in control theory. Its principal advantage is that it does not require a detailed system model.
 */
public class AlphaBetaFilter {
    private double position; // Estimated position
    private double velocity; // Estimated velocity
    private final double alpha; // Smoothing factor for position
    private final double beta; // Smoothing factor for velocity
    private final double dt; // Time step (delta time)

    // Constructor to initialize the filter
    public AlphaBetaFilter(double initialPosition, double initialVelocity, double alpha, double beta, double dt) {
        this.position = initialPosition;
        this.velocity = initialVelocity;
        this.alpha = alpha;
        this.beta = beta;
        this.dt = dt;
    }

    /**
     * Update the filter with a new measurement.
     *
     * @param measuredPosition The observed position.
     */
    public void update(double measuredPosition) {
        // Predict the next position
        double predictedPosition = position + velocity * dt;

        // Calculate the residual (error between measured and predicted positions)
        double residual = measuredPosition - predictedPosition;

        // Update position and velocity using alpha and beta gains
        position = predictedPosition + alpha * residual;
        velocity = velocity + beta * residual / dt;
    }

    /**
     * Get the current estimated position.
     *
     * @return Estimated position.
     */
    public double getPosition() {
        return position;
    }

    /**
     * Get the current estimated velocity.
     *
     * @return Estimated velocity.
     */
    public double getVelocity() {
        return velocity;
    }
}
