package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.extrapolation.ExtrapolationMethod;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;
import org.junit.Test;

import static org.junit.Assert.*;

public class CubicSplineTest {

    @Test
    public void testNaturalSplineInterpolation() {
        double[] x = new double[]{0.9, 1.3, 1.9, 2.1};
        double[] y = new double[]{1.3, 1.5, 1.85, 2.1};
        double[] xi = new double[]{1.15, 1.8, 2.0};
        CubicSpline cs = CubicSpline.newNaturalSpline(x, y);

        double yi = cs.evaluateAt(xi[0]);
        assertEquals(1.4307218309859153, yi, 1e-12);

        yi = cs.evaluateAt(xi[1]);
        assertEquals(1.755721830985916, yi, 1e-12);

        yi = cs.evaluateAt(xi[2]);
        assertEquals(1.9682218309859156, yi, 1e-12);

        double[] expected = {0.0, 21.438503269035536, 42.95370330964467, 64.62229689340103, 86.52098079187817,
                108.72645177664975, 131.31540661928935, 154.36454209137057, 177.95055496446702, 202.15014201015228,
                227.04, 252.68284377664975, 279.0854604670051, 306.2406552690355, 334.1412333807106, 362.78,
                392.15420158375633, 422.2788496243655, 453.1733968730964, 484.8572960812183, 517.35, 550.7171732791878,
                585.2093281624366, 621.1126457394247, 658.4708168934011, 697.0850423011844, 736.7459795871405,
                777.2442863756345, 818.3706202910322, 859.9156389576988, 901.67};

        x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
        cs = CubicSpline.newNaturalSpline(x, y);
        assertArrayEquals(expected, cs.evaluateAt(DoubleArrays.linSteps(0, 30)), 1e-12);

        assertEquals(392.15420158375633, cs.evaluateAt(16.0), 1e-12);
        assertEquals(29.746182686971242, cs.differentiate(16.0), 1e-12);
        assertEquals(29.746182686971242, cs.evaluateDerivativeAt(2, 16), 1e-12);
        assertEquals(377.4053741184433, cs.evaluateAntiDerivativeAt(2, 16.0), 1e-12);
        assertEquals(2947.3965035600677, cs.integrate(16.0), 1e-12);
        assertEquals(2947.3965035600677, cs.integrate(0.0, 16.0), 1e-12);
        assertEquals(1604.3556840203046, cs.integrate(11., 16.0), 1e-12);
    }

    @Test
    public void testParabolicallyTerminatedSplineInterpolation() {
        double[] x = new double[]{0.9, 1.3, 1.9, 2.1};
        double[] y = new double[]{1.3, 1.5, 1.85, 2.1};
        double[] xi = new double[]{1.15, 1.8, 2.0};
        CubicSpline cs = CubicSpline.newParabolicallyTerminatedSpline(x, y);

        double yi = cs.evaluateAt(xi[0]);
        assertEquals(1.4321022727272725, yi, 1e-12);

        yi = cs.evaluateAt(xi[1]);
        assertEquals(1.7594696969696972, yi, 1e-12);

        yi = cs.evaluateAt(xi[2]);
        assertEquals(1.9632575757575759, yi, 1e-12);

        double[] expected = {0.0, 20.134455464926596, 40.83992082653617, 62.11639608482872, 83.96388123980425,
                106.38237629146276, 129.37188123980425, 152.93239608482872, 177.06392082653616, 201.7664554649266,
                227.04, 252.89101146275146, 279.3517750081566, 306.4610328221859, 334.2575270908102, 362.78,
                392.0648895704187, 422.1394171615008, 453.0284999673736, 484.7570551821642, 517.35, 550.8445661337684,
                585.3272433713976, 620.8947775965198, 657.5965709624796, 695.4346818923328, 734.4091103860794,
                774.5198564437194, 815.7669200652529, 858.1503012506797, 901.67};

        x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
        cs = CubicSpline.newParabolicallyTerminatedSpline(x, y);
        assertArrayEquals(expected, cs.evaluateAt(DoubleArrays.linSteps(0, 30)), 1e-12);

        assertEquals(392.0648895704187, cs.evaluateAt(16.0), 1e-12);
        assertEquals(29.67555571506254, cs.differentiate(16.0), 1e-12);
        assertEquals(29.67555571506254, cs.evaluateDerivativeAt(2, 16), 1e-12);
        assertEquals(377.3576798332427, cs.evaluateAntiDerivativeAt(2, 16.0), 1e-12);
        assertEquals(2932.5665909262284, cs.integrate(16.0), 1e-12);
        assertEquals(2932.5665909262284, cs.integrate(0.0, 16.0), 1e-12);
        assertEquals(1605.0344416675732, cs.integrate(11., 16.0), 1e-12);
    }

    @Test
    public void testClampedSplineInterpolation() {
        double[] x = new double[]{0.9, 1.3, 1.9, 2.1};
        double[] y = new double[]{1.3, 1.5, 1.85, 2.1};
        double[] xi = new double[]{1.15, 1.8, 2.0};
        CubicSpline cs = CubicSpline.newClampedSpline(x, y, 2, 1);

        double yi = cs.evaluateAt(xi[0]);
        assertEquals(1.5100360576923078, yi, 1e-12);

        yi = cs.evaluateAt(xi[1]);
        assertEquals(1.7361111111111118, yi, 1e-12);

        yi = cs.evaluateAt(xi[2]);
        assertEquals(1.9814102564102565, yi, 1e-12);

        double[] expected = {0.0, 5.401105479452059, 17.012997260273984, 33.9485383561644, 55.32059178082194,
                80.24202054794523, 107.82568767123291, 137.1844561643836, 167.43118904109593, 197.67874931506856,
                227.04, 254.89402958904108, 281.6848284931507, 308.1226126027397, 334.91759780821917, 362.78,
                392.21514849315065, 422.90882630136986, 454.34192986301366, 485.99535561643836, 517.35,
                548.6513786301371, 583.2034838356165, 624.9080501268392, 673.8286410958905, 726.1906494165398,
                778.0525907661086, 825.4729808219179, 864.5103352612887, 891.2231697615423, 901.6700000000001};

        x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
        cs = CubicSpline.newClampedSpline(x, y, 2, 1);
        assertArrayEquals(expected, cs.evaluateAt(DoubleArrays.linSteps(0, 30)), 1e-12);

        assertEquals(392.21514849315065, cs.evaluateAt(16.0), 1e-12);
        assertEquals(30.15093041095891, cs.differentiate(16.0), 1e-12);
        assertEquals(30.15093041095891, cs.evaluateDerivativeAt(2, 16), 1e-12);
        assertEquals(377.37106748858446, cs.evaluateAntiDerivativeAt(2, 16.0), 1e-12);
        assertEquals(2765.246067488585, cs.integrate(16.0), 1e-12);
        assertEquals(2765.246067488585, cs.integrate(0.0, 16.0), 1e-12);
        assertEquals(1610.8140541552516, cs.integrate(11., 16.0), 1e-12);
    }

    @Test
    public void testNotAKnotSplineInterpolation() {
        double[] x = new double[]{0.9, 1.3, 1.9, 2.1};
        double[] y = new double[]{1.3, 1.5, 1.85, 2.1};
        double[] xi = new double[]{1.15, 1.8, 2.0};
        CubicSpline cs = CubicSpline.newNotAKnotSpline(x, y);

        double yi = cs.evaluateAt(xi[0]);
        assertEquals(1.4394531249999998, yi, 1e-12);

        yi = cs.evaluateAt(xi[1]);
        assertEquals(1.7593750000000004, yi, 1e-12);

        yi = cs.evaluateAt(xi[2]);
        assertEquals(1.9622916666666668, yi, 1e-12);

        double[] expected = {0.0, 20.514440000000018, 41.45489777777781, 62.844080000000034, 84.70469333333337,
                107.05944444444448, 129.93104000000002, 153.3421866666667, 177.31559111111113, 201.87396, 227.04,
                252.83641777777777, 279.28592, 306.4112133333333, 334.2350044444444, 362.78, 392.0707644444444,
                422.13929333333334, 453.01944000000003, 484.7450577777778, 517.35, 550.87058, 585.3529511111111,
                620.8457266666667, 657.39752, 695.0569444444445, 733.8726133333333, 773.89314, 815.1671377777777,
                857.7432200000001, 901.67};

        x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
        cs = CubicSpline.newNotAKnotSpline(x, y);
        assertArrayEquals(expected, cs.evaluateAt(DoubleArrays.linSteps(0, 30)), 1e-12);

        assertEquals(392.0707644444444, cs.evaluateAt(16.0), 1e-12);
        assertEquals(29.674004444444453, cs.differentiate(16.0), 1e-12);
        assertEquals(29.674004444444453, cs.evaluateDerivativeAt(2, 16), 1e-12);
        assertEquals(377.36197907407404, cs.evaluateAntiDerivativeAt(2, 16.0), 1e-12);
        assertEquals(2936.883854074074, cs.integrate(16.0), 1e-12);
        assertEquals(2936.883854074074, cs.integrate(0.0, 16.0), 1e-12);
        assertEquals(1604.869493148148, cs.integrate(11., 16.0), 1e-12);
    }

    @Test
    public void testAkimaSplineInterpolation() {
        double[] x = {0.5, 0.9, 1.3, 1.9, 2.1, 2.2};
        double[] y = {1.0, 1.3, 1.5, 1.85, 2.1, 2.4};
        double[] xi = new double[]{1.15, 1.8, 2.0};
        CubicSpline cs = CubicSpline.newAkimaSpline(x, y);

        double yi = cs.evaluateAt(xi[0]);
        //1.425865589488636
        assertEquals(1.4280628551136363, yi, 1e-12);

        yi = cs.evaluateAt(xi[1]);
        assertEquals(1.788720538720539, yi, 1e-12);

        yi = cs.evaluateAt(xi[2]);
        assertEquals(1.9299242424242427, yi, 1e-12);

        double[] expected = {0.0, 20.0376, 40.66773333333334, 61.8904, 83.7056, 106.11333333333333, 129.11359999999996,
                152.70639999999997, 176.8917333333333, 201.66959999999997, 227.04, 252.89567093511613,
                279.2366128053484, 306.2237192080226, 334.0178837404645, 362.78, 392.4458809357619, 422.80077482065684,
                453.7797282376708, 485.31778776978985, 517.35, 550.5467494784905, 585.2404776618868, 620.9308, 657.6892,
                695.5633333333334, 734.5532000000001, 774.6587999999999, 815.8801333333333, 858.2171999999999, 901.67};

        x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
        cs = CubicSpline.newAkimaSpline(x, y);
        assertArrayEquals(expected, cs.evaluateAt(DoubleArrays.linSteps(0, 30)), 1e-12);

        assertEquals(392.4458809357619, cs.evaluateAt(16.0), 1e-12);
        assertEquals(30.021212979830768, cs.differentiate(16.0), 1e-12);
        assertEquals(30.021212979830768, cs.evaluateDerivativeAt(2, 16), 1e-12);
        assertEquals(377.5528163297443, cs.evaluateAntiDerivativeAt(2, 16.0), 1e-12);
        assertEquals(2930.356170519625, cs.integrate(16.0), 1e-12);
        assertEquals(2930.356170519625, cs.integrate(0.0, 16.0), 1e-12);
        assertEquals(1604.5998481745492, cs.integrate(11., 16.0), 1e-12);
    }

    @Test
    public void testGetSegmentNotAKnotSplineInterpolation() {
        double[] x = new double[]{0.9, 1.3, 1.9, 2.1};
        double[] y = new double[]{1.3, 1.5, 1.85, 2.1};
        double[] xi = new double[]{1.15, 1.8, 2.0};
        CubicSpline cs = CubicSpline.newNotAKnotSpline(x, y);

        double yi = cs.getFirstSegment().evaluateAt(xi[0]);
        assertEquals(1.4394531249999998, yi, 1e-12);

        yi = cs.getSegment(1).evaluateAt(xi[1]);
        assertEquals(1.7593750000000004, yi, 1e-12);

        yi = cs.getLastSegment().evaluateAt(xi[2]);
        assertEquals(1.9622916666666668, yi, 1e-12);
    }

    @Test
    public void testNaturalSplineExtrapolateLeft() {
        double[] x = new double[]{0.9, 1.3, 1.9, 2.1};
        double[] y = new double[]{1.3, 1.5, 1.85, 2.1};
        CubicSpline cs = CubicSpline.newNaturalSpline(x, y);
        double left = -0.5;

        cs.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_END_POINT);
        double yi = cs.evaluateAt(left);
        assertEquals("Natural Spline ClampToEndPoint lower bound extrapolation", 1.3, yi, 0.0);

        cs.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_NAN);
        yi = cs.evaluateAt(left);
        assertTrue("Natural Spline ClampToEndNaN lower bound extrapolation", Double.isNaN(yi));

        cs.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_ZERO);
        yi = cs.evaluateAt(left);
        assertEquals("Natural Spline ClampToZero lower bound extrapolation", 0.0, yi, 0.0);

        cs.setExtrapolationMethod(ExtrapolationMethod.LINEAR);
        yi = cs.evaluateAt(left);
        assertEquals("Natural Spline Linear lower bound extrapolation", 0.5474178403755874, yi,
                1e-12);

        cs.setExtrapolationMethod(ExtrapolationMethod.NATURAL);
        yi = cs.evaluateAt(left);
        assertEquals("Natural Spline Natural lower bound extrapolation", 1.1915492957746414, yi,
                1e-12);

        cs.setExtrapolationMethod(ExtrapolationMethod.THROW);
        IndexOutOfBoundsException exception = assertThrows(IndexOutOfBoundsException.class, () -> {
            cs.evaluateAt(left);
        });
        assertEquals("x = -0.5000 is smaller than every number in x[]", exception.getMessage());
    }

    @Test
    public void testNaturalSplineExtrapolateRight() {
        double[] x = new double[]{0.9, 1.3, 1.9, 2.1};
        double[] y = new double[]{1.3, 1.5, 1.85, 2.1};
        CubicSpline cs = CubicSpline.newNaturalSpline(x, y);
        double right = 3.0;

        cs.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_END_POINT);
        double yi = cs.evaluateAt(right);
        assertEquals("Natural Spline ClampToEndPoint upper bound extrapolation", 2.1, yi, 0.0);

        cs.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_NAN);
        yi = cs.evaluateAt(right);
        assertTrue("Natural Spline ClampToEndNaN upper bound extrapolation", Double.isNaN(yi));

        cs.setExtrapolationMethod(ExtrapolationMethod.CLAMP_TO_ZERO);
        yi = cs.evaluateAt(right);
        assertEquals("Natural Spline ClampToZero upper bound extrapolation", 0.0, yi, 0.0);

        cs.setExtrapolationMethod(ExtrapolationMethod.LINEAR);
        yi = cs.evaluateAt(right);
        assertEquals("Natural Spline Linear upper bound extrapolation", 3.306338028169013, yi,
                1e-12);

        cs.setExtrapolationMethod(ExtrapolationMethod.NATURAL);
        yi = cs.evaluateAt(right);
        assertEquals("Natural Spline Natural upper bound extrapolation", 1.6592429577464949, yi,
                1e-12);

        cs.setExtrapolationMethod(ExtrapolationMethod.THROW);
        IndexOutOfBoundsException exception = assertThrows(IndexOutOfBoundsException.class, () -> {
            cs.evaluateAt(right);
        });
        assertEquals("x = 3.0000 is bigger than every number in x[]", exception.getMessage());
    }
}
