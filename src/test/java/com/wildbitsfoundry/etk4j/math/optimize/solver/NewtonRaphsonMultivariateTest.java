package com.wildbitsfoundry.etk4j.math.optimize.solver;

import com.wildbitsfoundry.etk4j.math.function.MultivariateFunction;
import com.wildbitsfoundry.etk4j.math.optimize.OptimizerStatusType;
import org.junit.Test;

import static org.junit.Assert.*;

public class NewtonRaphsonMultivariateTest {

    @Test
    public void testNewtonRaphsonMethodPreComputedJacobianDefaultConditions() {
        // Define the system of equations
        MultivariateFunction[] functions = {
                x -> x[0] + x[1] - 3 * x[2] + x[3] - 2,
                x -> -5 * x[0] + 3 * x[1] - 4 * x[2] + x[3],
                x -> x[0] + 2 * x[2] - x[3] - 1,
                x -> x[0] + 2 * x[1] - 12
        };
        // Define Jacobian
        MultivariateFunction[][] jacobian = {
                {x -> 1, x -> 1, x -> -3, x -> 1},
                {x -> -5, x -> 3, x -> -4, x -> 1},
                {x -> 1, x -> 0, x -> 2, x -> -1},
                {x -> 1, x -> 2, x -> 0, x -> 0}
        };
        //Initial guess
        double[] x0 = {1, 5, 5, 10};

        // Solve using the Newton-Raphson method
        SolverResults<double[]> nr = new NewtonRaphsonMultivariate(functions, x0)
                .jacobian(jacobian)
                .solve();
        double[] expected = {1.2941176470588236, 5.352941176470588, 4.9411764705882355, 10.176470588235293};
        assertArrayEquals(expected, nr.getValue(), 1e-12);
        assertEquals(1, nr.getNumberOfIterations());
        assertEquals("Converged", nr.getSolverStatus());
        assertEquals(OptimizerStatusType.CONVERGED, nr.getOptimizerStatusType());
        assertEquals(3.076740298213702E-15, nr.getError(), 1e-12);
        assertTrue(nr.hasConverged());
    }

    @Test
    public void testNewtonRaphsonMethodPreComputedJacobian() {
        // Define the system of equations
        MultivariateFunction[] functions = {
                x -> x[0] + x[1] - 3 * x[2] + x[3] - 2,
                x -> -5 * x[0] + 3 * x[1] - 4 * x[2] + x[3],
                x -> x[0] + 2 * x[2] - x[3] - 1,
                x -> x[0] + 2 * x[1] - 12
        };
        // Define Jacobian
        MultivariateFunction[][] jacobian = {
                {x -> 1, x -> 1, x -> -3, x -> 1},
                {x -> -5, x -> 3, x -> -4, x -> 1},
                {x -> 1, x -> 0, x -> 2, x -> -1},
                {x -> 1, x -> 2, x -> 0, x -> 0}
        };
        //Initial guess
        double[] x0 = {1, 5, 5, 10};

        // Solve using the Newton-Raphson method
        SolverResults<double[]> nr = new NewtonRaphsonMultivariate(functions, x0)
                .jacobian(jacobian)
                .tolerance(1e-6)
                .iterationLimit(100)
                .differentiationStepSize(1e-6)
                .lineSearchArmijoParameter(1e-4)
                .lineSearchInitialStepSize(1)
                .lineSearchStepSizeReductionFactor(0.9)
                .lineSearchIterationLimit(100)
                .solve();
        double[] expected = {1.2941176470588236, 5.352941176470588, 4.9411764705882355, 10.176470588235293};
        assertArrayEquals(expected, nr.getValue(), 1e-12);
        assertEquals(1, nr.getNumberOfIterations());
        assertEquals("Converged", nr.getSolverStatus());
        assertEquals(OptimizerStatusType.CONVERGED, nr.getOptimizerStatusType());
        assertEquals(3.076740298213702E-15, nr.getError(), 1e-12);
        assertTrue(nr.hasConverged());
    }

    @Test
    public void testNewtonRaphsonMethodApproximatedJacobian() {
        // Define the system of equations
        MultivariateFunction[] functions = {
                x -> x[0] + x[1] - 3 * x[2] + x[3] - 2,
                x -> -5 * x[0] + 3 * x[1] - 4 * x[2] + x[3],
                x -> x[0] + 2 * x[2] - x[3] - 1,
                x -> x[0] + 2 * x[1] - 12
        };
        //Initial guess
        double[] x0 = {1, 5, 5, 10};
        // Solve using the Newton-Raphson method

        SolverResults<double[]> nr = new NewtonRaphsonMultivariate(functions, x0)
                .jacobian(null)
                .tolerance(1e-6)
                .iterationLimit(100)
                .differentiationStepSize(1e-6)
                .lineSearchArmijoParameter(1e-4)
                .lineSearchInitialStepSize(1)
                .lineSearchStepSizeReductionFactor(0.9)
                .lineSearchIterationLimit(100)
                .solve();
        double[] expected = {1.29411764689683, 5.352941176655849, 4.941176470937592, 10.176470588632691};
        assertArrayEquals(expected, nr.getValue(), 1e-12);
        assertEquals(1, nr.getNumberOfIterations());
        assertEquals("Converged", nr.getSolverStatus());
        assertEquals(OptimizerStatusType.CONVERGED, nr.getOptimizerStatusType());
        assertEquals(7.682979784690616E-10, nr.getError(), 1e-12);
        assertTrue(nr.hasConverged());
    }

    @Test
    public void testNewtonRaphsonMethodPreComputedJacobianDefaultConditionsCosineAndSine() {
        // Define the system of equations
        MultivariateFunction[] functions = {
                x -> 3 * x[0] - Math.cos(x[1] * x[2]) - 3.0/2.0,
                x -> 4 * x[0] * x[0] -625 * x[1] * x[1] + 2 * x[2] - 1,
                x -> 20 * x[2] + Math.exp(-x[0] * x[1]) + 9
        };
        // Define Jacobian
        MultivariateFunction[][] jacobian = {
                {x -> 3, x -> x[2] * Math.sin(x[1] * x[2]), x -> x[1] * Math.sin(x[1] * x[2])},
                {x -> 8 * x[0], x -> -1250 * x[1], x -> 2},
                {x -> -x[1] * Math.exp(-x[0] * x[1]), x -> -x[0] * Math.exp(-x[0] * x[1]), x -> 20}
        };
        //Initial guess
        double[] x0 = {1, 1, 1};

        // Solve using the Newton-Raphson method
        SolverResults<double[]> nr = new NewtonRaphsonMultivariate(functions, x0)
                .jacobian(jacobian)
                .solve();
        double[] expected = {0.8332816138167236, 0.035334617159854045, -0.49854927776854513};
        assertArrayEquals(expected, nr.getValue(), 1e-12);
        assertEquals(8, nr.getNumberOfIterations());
        assertEquals("Converged", nr.getSolverStatus());
        assertEquals(OptimizerStatusType.CONVERGED, nr.getOptimizerStatusType());
        assertEquals(4.498299174920495E-8, nr.getError(), 1e-12);
        assertTrue(nr.hasConverged());
    }
}
