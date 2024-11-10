package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;

public class FiltersTest {

    @Test
    public void testLptolp() {
        double[] num = {1};
        double[] den = {3, 2, 1};

        TransferFunction tf = Filters.lpToLp(num, den, 10.0);
        assertArrayEquals(new double[]{33.33333333333333}, tf.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 6.666666666666666, 33.33333333333333},
                tf.getDenominator().getCoefficients(), 1e-12);
    }

    @Test
    public void testLpTohp() {
        double[] num = {1};
        double[] den = {3, 2, 1};

        TransferFunction tf = Filters.lpToHp(num, den, 10.0);
        assertArrayEquals(new double[]{1.0, 0.0, 0.0}, tf.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 20.0, 300.0},
                tf.getDenominator().getCoefficients(), 1e-12);
    }

    @Test
    public void testLpTobp() {
        double[] num = {1};
        double[] den = {1, 2, 1};

        TransferFunction tf = Filters.lpToBp(num, den, 10.0, 100.0);
        assertArrayEquals(new double[]{10000.0, 0.0, 0.0}, tf.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 200.0, 10200.0, 20000.0, 10000.0},
                tf.getDenominator().getCoefficients(), 1e-12);
    }

    @Test
    public void testLpTobs() {
        double[] num = {1};
        double[] den = {1, 2, 1};

        TransferFunction tf = Filters.lpToBs(num, den, 10.0, 100.0);
        assertArrayEquals(new double[]{1.0, 0.0, 200.0, 0.0, 10000.0}, tf.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 200.0, 10200.0, 20000.0, 10000.0},
                tf.getDenominator().getCoefficients(), 1e-12);
    }
}
