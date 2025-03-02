package com.wildbitsfoundry.etk4j.signal.stoz;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.signal.filters.ButterWorth;
import org.junit.Test;

import static com.wildbitsfoundry.etk4j.signal.stoz.BilinearTransform.transform;
import static org.junit.Assert.assertArrayEquals;

public class BilinearTransformTest {

    @Test
    public void testTransformBandpass() {
        TransferFunction filts = ButterWorth.newBandpass(4, 2 * Math.PI * 7, 2 * Math.PI * 13);
        double[][] filtz = transform(filts.getNumeratorCoefficients(), filts.getDenominatorCoefficients(), 100);

        double[] num = {5.705645409457357E-4, 0.0, -0.002282258163782943, 3.5309398223506414E-19, 0.003423387245674412,
                -3.5309398223506414E-19, -0.002282258163782943, 0.0, 5.705645409457357E-4};

        double[] den = {0.9999999999999999, -5.935553957855101, 16.373268483011213, -27.189994073125547,
                29.643566690552934, -21.702719846722392, 10.431386855712248, -3.0188859533519223, 0.4064602384544566};

        assertArrayEquals(num, filtz[0], 1e-12);
        assertArrayEquals(den, filtz[1], 1e-12);
    }
}
