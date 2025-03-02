package com.wildbitsfoundry.etk4j.math.specialfunction;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class GammaTest {

    @Test
    public void testGammaFunction() {
        double[] odd = {0};
        double[] even = {0};
        double solution = Gamma.gamma(0.1, odd, even);
        assertEquals(9.51350769866873, solution, 1e-12);
    }

    @Test
    public void testRecipGammaFunction() {
        double[] odd = {0};
        double[] even = {0};
        double solution = Gamma.recipGamma(0.1, odd, even);
        assertEquals(0.9357787209128726, solution, 1e-12);
    }

    @Test
    public void testLogGammaFunction() {
        double solution = Gamma.logGamma(0.1);
        assertEquals(2.252712651734206, solution, 1e-12);
    }
}
