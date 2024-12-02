package com.wildbitsfoundry.etk4j.math.calculus.odesolvers;

import com.wildbitsfoundry.etk4j.math.functions.BivariateFunction;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertArrayEquals;

public class RungeKuttaTest {

    @Test
    public void testRungeKutta23() {
        BivariateFunction func = (t, x) -> -x;
        RungeKutta rungeKutta = new RungeKutta23(func, 0.0, 1.0, 10.0, 1.0, 0.001,
                Math.exp(-6), null);

        List<Double> tValues = new ArrayList<>();
        List<Double> yValues = new ArrayList<>();
        while (!rungeKutta.status.equals("finished")) {
            rungeKutta.step();
            tValues.add(rungeKutta.t);
            yValues.add(rungeKutta.y[0]);
        }

        double[] tSolution = {0.0326443353286873, 0.35908768861556034, 0.9288876265576596, 1.6469479545357713,
                2.6156116724039244, 3.6156116724039244, 4.615611672403924, 5.615611672403924, 6.615611672403924,
                7.615611672403924, 8.615611672403924, 9.615611672403924, 10.0};

        double[] ySolution = {0.9678826930655442, 0.6978834512680837, 0.392003252913781, 0.1873926074634615,
                0.06540125734101757, 0.02180041911367253, 0.007266806371224182, 0.002422268790408061,
                8.074229301360204E-4, 2.691409767120068E-4, 8.971365890400227E-5, 2.9904552968000765E-5,
                2.033578448145516E-5};

        assertArrayEquals(tSolution, tValues.stream().mapToDouble(Double::doubleValue).toArray(), 1e-12);
        assertArrayEquals(ySolution, yValues.stream().mapToDouble(Double::doubleValue).toArray(), 1e-12);
    }

    @Test
    public void testRungeKutta23MultipleInitialConditions() {
        ODESystemOfEquations odeSystemOfEquations = (t, y) -> {
            double dxdt = y[0] - y[1];
            double dydt = y[0] + y[1];
            return new double[] {dxdt, dydt};
        };
        RungeKutta rungeKutta = new RungeKutta23(odeSystemOfEquations, 0.0, new double[] {1, 0}, 10.0, 1.0, 0.001,
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

        double[] tSolution = {0.025976025597427347, 0.28573628157170083, 0.6423974632107396, 0.9415602755161425,
                1.2407230878215454, 1.517279443855995, 1.7715083511313474, 2.021787100256841, 2.29097421924673,
                2.5514070945230576, 2.7905309430240197, 3.0126104081129306, 3.1857583889400614, 3.358906369767192,
                3.563997453026649, 3.8057525582706306, 4.05090314219944, 4.277521220595271, 4.489607069883245,
                4.6577675037614314, 4.825927937639618, 4.9780007440115925, 5.182210064824636, 5.421904791058065,
                5.660287066661956, 5.88084070681308, 6.086840762766752, 6.245191922304215, 6.403543081841679,
                6.552534026960952, 6.7548446113692995, 6.993402989396838, 7.230866785247203, 7.450639053591104,
                7.655946482400699, 7.813594830473213, 7.971243178545727, 8.117344093508647, 8.317274652391955,
                8.554840682474776, 8.792688011927234, 9.0128248507089, 9.218633356666684, 9.377736872424684,
                9.52089620353776, 9.65618216873147, 9.846268246481957, 10.0};

        double[] ySolution0 = {1.0259701831225174, 1.2776035178859115, 1.5267550915211654, 1.5163127433229293,
                1.129856646197454, 0.25140486100338855, -1.1706776106990222, -3.3020115943796675, -6.552776469076043,
                -10.725013995834031, -15.41073920810532, -20.33853485998865, -24.373766691221242, -28.34124295452492,
                -32.52272281611676, -35.79963073333264, -35.74685313257538, -30.82683129280548, -20.100592247483277,
                -6.0696363053903575, 13.98169102504604, 38.18775616304511, 81.13082982480888, 148.7128818678593,
                235.66673642208997, 333.23109350974744, 436.8143851935606, 521.6862994454893, 607.532293124006,
                685.0106992865498, 775.6835225602056, 839.344527619152, 820.8462522725175, 690.0179629245631,
                428.2482878009952, 109.03692017439732, -334.9230218116534, -875.5305695946774, -1845.669436294149,
                -3387.276696489381, -5385.505416681393, -7632.654713679352, -10026.323060984238, -12004.132180751249,
                -13810.64538648963, -15469.332061476902, -17554.20169595952, -18838.970153323444};

        double[] ySolution1 = {0.02665662197817555, 0.37515343355827274, 1.1407960593646902, 2.078914915154872,
                3.2851597315719747, 4.577379565118072, 5.797551084526619, 6.845814193180313, 7.494510471488281,
                7.212614936750099, 5.6777033023639545, 2.679156690185075, -1.0281111935324643, -6.197521909276459,
                -14.536889066248099, -27.899308603442968, -45.70545510166953, -65.96119404243086, -87.7634441762113,
                -106.36305972168594, -125.28245468783436, -141.72167429578568, -160.76108836518648, -174.0763962466806,
                -170.21429656456093, -142.86326107330848, -88.1191636942311, -21.254931925323348, 71.78889956029016,
                187.07842083163416, 392.92284921345833, 717.5800846971679, 1135.1666470317791, 1603.578622691247,
                2100.920783997546, 2508.0988359639473, 2920.2626567826746, 3287.450183903881, 3723.5799897544484,
                4040.652354197807, 3971.1180706908385, 3365.5047092084988, 2132.032837865074, 605.8951467207048,
                -1284.2212340228366, -3589.4499602024716, -7798.406381055847, -12122.588618481495};

        assertArrayEquals(tSolution, tValues.stream().mapToDouble(Double::doubleValue).toArray(), 1e-12);
        assertArrayEquals(ySolution0, yValues0.stream().mapToDouble(Double::doubleValue).toArray(), 1e-12);
        assertArrayEquals(ySolution1, yValues1.stream().mapToDouble(Double::doubleValue).toArray(), 1e-12);
    }

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
