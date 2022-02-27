package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class ZeroPoleGainTest {
    @Test
    public void testEvaluateAt() {
        double[] numerator = {1, 2, 1};
        double[] denominator = {1, 2, 2, 1};
        Complex[] zeros = new Polynomial(numerator).calculateRoots();
        Complex[] poles = new Polynomial(denominator).calculateRoots();
        double gain = 2.0;
        ZeroPoleGain zpk = new ZeroPoleGain(zeros, poles, gain);
        TransferFunction tf = new TransferFunction(numerator, denominator).multiply(gain);
        assertEquals(tf.evaluateAt(100), zpk.evaluateAt(100));
    }
}
