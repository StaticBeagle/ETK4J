package com.wildbitsfoundry.etk4j.math.calculus.odesolvers;

import com.wildbitsfoundry.etk4j.math.functions.BivariateFunction;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertArrayEquals;

public class RungeKuttaTest {

    @Test
    public void testRungeKutta45() {
        BivariateFunction func = (t, x) -> -x;
        RungeKutta rungeKutta = new RungeKutta45(func, 0.0, 1.0, 10.0, 1.0, 0.001,
                Math.exp(-6), null);

        List<Double> tValues = new ArrayList<>();
        List<Double> yValues = new ArrayList<>();
        while (!rungeKutta.status.equals("finished")) {
            rungeKutta.step();
            tValues.add(rungeKutta.t);
            yValues.add(rungeKutta.y[0]);
        }

        double[] tSolution = {0.1283171479634216, 1.1283171479634215, 2.1283171479634215, 3.1283171479634215,
                4.128317147963422, 5.128317147963422, 6.128317147963422, 7.128317147963422, 8.128317147963422,
                9.128317147963422, 10.0};

        double[] ySolution = {0.879574381033538, 0.3239765636806864, 0.11933136762238628, 0.043953720407578944,
                0.01618962035012491, 0.005963176828962677, 0.002196436798667919, 8.090208875093502E-4,
                2.9798936023261037E-4, 1.097594143523445E-4, 4.5927433621121034E-5};

        assertArrayEquals(tSolution, tValues.stream().mapToDouble(Double::doubleValue).toArray(), 1e-12);
        assertArrayEquals(ySolution, yValues.stream().mapToDouble(Double::doubleValue).toArray(), 1e-12);
    }

    @Test
    public void testRungeKutta45MultipleInitialConditions() {
        ODESystemOfEquations odeSystemOfEquations = (t, y) -> {
            double dxdt = y[0] - y[1];
            double dydt = y[0] + y[1];
            return new double[] {dxdt, dydt};
        };
        RungeKutta rungeKutta = new RungeKutta45(odeSystemOfEquations, 0.0, new double[] {1, 0}, 10.0, 1.0, 0.001,
                Math.exp(-6), null);

        List<Double> tValues = new ArrayList<>();
        List<Double> yValues0 = new ArrayList<>();
        List<Double> yValues1 = new ArrayList<>();
        while (!rungeKutta.status.equals("finished")) {
            rungeKutta.step();
            tValues.add(rungeKutta.t);
            yValues0.add(rungeKutta.y[0]);
            yValues1.add(rungeKutta.y[1]);
        }

        double[] tSolution = {0.11187772386872015, 1.007930843193444, 1.999271727606395, 2.8448211368811007,
                3.684921261666479, 4.568816163110601, 5.382471317367101, 6.232827184027202, 7.049004403889039,
                7.886515731920654, 8.711704819712274, 9.539576468671115, 10.0};

        double[] ySolution0 = {1.1113842514617591, 1.4610914926227114, -3.069122887043963, -16.429902353989913,
                -34.035760621658426, -13.677393049997438, 135.0117027757588, 507.3684317131492, 827.0322316967461,
                -89.79023363872113, -4584.94515450421, -13758.728742928337, -18401.94942485325};

        double[] ySolution1 = {0.12486051574905874, 2.3166396846717032, 6.70775547689553, 5.011061651185726,
                -20.594892705054463, -95.27367878475552, -170.01846262271692, -24.983562453576752, 797.212040125511,
                2651.9899373783214, 3954.3666464557764, -1606.8576979577356, -11969.226607049524};

        assertArrayEquals(tSolution, tValues.stream().mapToDouble(Double::doubleValue).toArray(), 1e-12);
        assertArrayEquals(ySolution0, yValues0.stream().mapToDouble(Double::doubleValue).toArray(), 1e-12);
        assertArrayEquals(ySolution1, yValues1.stream().mapToDouble(Double::doubleValue).toArray(), 1e-12);
    }
}
