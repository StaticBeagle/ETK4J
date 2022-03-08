package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.util.DoubleArrays;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;

public class NearestNeighborTest {

    @Test
    public void testInterpolation() {
        double[] x = DoubleArrays.linSteps(0, 4, 0.5);
        double[] y = Arrays.stream(x).map(v -> v * v).toArray();

        NearestNeighbor nh = NearestNeighbor.newNearestNeighbor(x, y);
        double[] xi = DoubleArrays.linSteps(0, 4, 0.05);
        double[] expected = {0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1.0,
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25,
                4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 6.25, 6.25, 6.25, 6.25, 6.25, 6.25, 6.25, 6.25, 6.25,
                6.25, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 12.25, 12.25, 12.25, 12.25, 12.25, 12.25, 12.25,
                12.25, 12.25, 12.25, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0};
        assertArrayEquals(expected, nh.evaluateAt(xi), 1e-12);
    }
}
