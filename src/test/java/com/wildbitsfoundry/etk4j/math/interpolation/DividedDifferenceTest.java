package com.wildbitsfoundry.etk4j.math.interpolation;

import org.junit.Test;

import static com.wildbitsfoundry.etk4j.math.interpolation.DividedDifference.computeDividedDifferenceTable;
import static com.wildbitsfoundry.etk4j.math.interpolation.DividedDifference.newtonPolynomial;
import static org.junit.Assert.assertEquals;

public class DividedDifferenceTest {

    @Test
    public void testInterpolation() {
        // Given data points
        double[] x = {8, 9, 9.5, 11};
        double[] y = {2.079442, 2.197225, 2.251292, 2.397895};

        DividedDifference dividedDifference = new DividedDifference(x, y);

        double result = dividedDifference.evaluateAt(9.2);
        assertEquals(2.21920816, result, 1e-12);
    }

    @Test
    public void testInterpolationUsingStaticMethods() {
        // Given data points
        double[] x = {8, 9, 9.5, 11};
        double[] y = {2.079442, 2.197225, 2.251292, 2.397895};

        // Compute the divided difference table
        double[][] table = computeDividedDifferenceTable(x, y);

        // Evaluate the Newton interpolating polynomial at a specific value
        double value = 9.2;
        double result = newtonPolynomial(x, table, value);
        assertEquals(2.21920816, result, 1e-12);
    }
}
