package com.wildbitsfoundry.etk4j.math.laplace;

import com.wildbitsfoundry.etk4j.math.functions.ComplexUnivariateFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class LaplaceTest {

    @Test
    public void testLaplace() {
        UnivariateFunction timeDomain = t -> Math.exp(-2.0 * t);
        UnivariateFunction frequencyDomain = s -> 1.0 / (s + 2);

        double frequency = 0;
        assertEquals(frequencyDomain.evaluateAt(frequency), Laplace.transform(timeDomain, frequency), 1e-3);

        frequency = 1;
        assertEquals(frequencyDomain.evaluateAt(frequency), Laplace.transform(timeDomain, frequency), 1e-3);

        frequency = 2;
        assertEquals(frequencyDomain.evaluateAt(frequency), Laplace.transform(timeDomain, frequency), 1e-3);

        frequency = 8;
        assertEquals(frequencyDomain.evaluateAt(frequency), Laplace.transform(timeDomain, frequency), 1e-3);
    }

    @Test
    public void testTalbot() {
        InverseLaplaceTransformTalbot lp = new InverseLaplaceTransformTalbot();

        UnivariateFunction timeDomain = t -> Math.exp(-2.0 * t);
        ComplexUnivariateFunction frequencyDomain = s -> s.add(2.0).invert();

        double time = 0.1;
        assertEquals(timeDomain.evaluateAt(time), lp.inverseTransform(frequencyDomain, time), 1e-3);

        time = 0.2;
        assertEquals(timeDomain.evaluateAt(time), lp.inverseTransform(frequencyDomain, time), 1e-3);

        time = 0.5;
        assertEquals(timeDomain.evaluateAt(time), lp.inverseTransform(frequencyDomain, time), 1e-3);

        time = 1.0;
        assertEquals(timeDomain.evaluateAt(time), lp.inverseTransform(frequencyDomain, time), 1e-3);

        time = 10.0;
        assertEquals(timeDomain.evaluateAt(time), lp.inverseTransform(frequencyDomain, time), 1e-3);

        time = 100.0;
        assertEquals(timeDomain.evaluateAt(time), lp.inverseTransform(frequencyDomain, time), 1e-3);
    }

    @Test
    public void testStehfest() {
        InverseLaplaceTransformStehfest lp = new InverseLaplaceTransformStehfest();

        UnivariateFunction timeDomain = t -> Math.exp(-2.0 * t);
        UnivariateFunction frequencyDomain = s -> 1.0 / (s + 2);

        double time = 0.1;
        assertEquals(timeDomain.evaluateAt(time), lp.inverseTransform(frequencyDomain, time), 1e-3);

        time = 0.2;
        assertEquals(timeDomain.evaluateAt(time), lp.inverseTransform(frequencyDomain, time), 1e-3);

        time = 0.5;
        assertEquals(timeDomain.evaluateAt(time), lp.inverseTransform(frequencyDomain, time), 1e-3);

        time = 1.0;
        assertEquals(timeDomain.evaluateAt(time), lp.inverseTransform(frequencyDomain, time), 1e-3);

        time = 100.0;
        assertEquals(timeDomain.evaluateAt(time), lp.inverseTransform(frequencyDomain, time), 1e-3);
    }
}
