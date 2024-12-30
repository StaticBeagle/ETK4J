package com.wildbitsfoundry.etk4j.signals.smoothing;

import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;

public class AlphaBetaFilterTest {

    @Test
    public void testUpdate() {
        // Example usage of the Alpha-Beta Filter
        double initialPosition = 0.0;
        double initialVelocity = 0.0;
        double alpha = 0.85; // Position smoothing factor
        double beta = 0.005; // Velocity smoothing factor
        double dt = 1.0; // Time step (e.g., 1 second)

        // Initialize the filter
        AlphaBetaFilter filter = new AlphaBetaFilter(initialPosition, initialVelocity, alpha, beta, dt);

        // Example measurements
        double[] measurements = {0.0, 1.2, 2.8, 4.1, 5.9, 7.4, 9.2};

        double[] expectedPosition = {0.0, 1.02, 2.5339, 3.8673154999999997, 5.5984912475, 7.1346751533875, 8.896429363503687};
        double[] expectedVelocity = {0.0, 0.006, 0.01487, 0.022626149999999998, 0.03267644175, 0.04152060330375, 0.05163962452029375};
        double[] actualPosition = new double[measurements.length];
        double[] actualVelocity = new double[measurements.length];

        for(int i = 0; i < measurements.length; i++) {
            filter.update(measurements[i]);
            actualPosition[i] = filter.getPosition();
            actualVelocity[i] = filter.getVelocity();
        }
        assertArrayEquals(expectedPosition, actualPosition, 1e-12);
        assertArrayEquals(expectedVelocity, actualVelocity, 1e-12);
    }
}
