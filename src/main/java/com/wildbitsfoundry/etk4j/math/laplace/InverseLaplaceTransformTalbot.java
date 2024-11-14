package com.wildbitsfoundry.etk4j.math.laplace;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.functions.ComplexUnivariateFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

/**
 * The {@code InverseLaplaceTransformTalbot} implements the inverse laplace transform using Talbot's method.
 * <pre>
 *     References:
 *     Tucker McClure (2022). Numerical Inverse Laplace Transform
 *     (https://www.mathworks.com/matlabcentral/fileexchange/39035-numerical-inverse-laplace-transform),
 *     MATLAB Central File Exchange. Retrieved February 16, 2022.
 * </pre>
 * See <a href="http://www.columbia.edu/~ww2040/UnifiedDraft.pdf">Unified Framework pages 17-18</a>
 */
public class InverseLaplaceTransformTalbot {

    private int M;
    private Complex[] delta;
    private Complex[] gamma;

    public InverseLaplaceTransformTalbot() {
        this(64);
    }

    /**
     * Construct an Inverse Laplace Transform using Talbot's method.
     * @param M The number of times to sum each time sample used in
     * {@link #inverseTransform(ComplexUnivariateFunction, double)}
     */
    public InverseLaplaceTransformTalbot(int M) {
        this.M = M;
        delta = ComplexArrays.zeros(M);
        delta[0] = Complex.fromReal(2.0 * M / 5.0);
        for (int k = 1; k < M; ++k) {
            delta[k] = new Complex(1.0 / Math.tan(k * Math.PI / M), 1.0);
            delta[k].multiplyEquals(2.0 * k * Math.PI / 5.0);
        }
        gamma = ComplexArrays.zeros(M);
        gamma[0] = delta[0].exp().multiply(0.5);
        for (int k = 1; k < M; ++k) {
            double cotTerm = 1.0 / Math.tan(k * Math.PI / M);
            gamma[k] = new Complex(1.0, (k * Math.PI / M) * (1.0 + cotTerm * cotTerm));
            gamma[k].subtractEquals(Complex.fromImaginary(cotTerm));
            gamma[k].multiplyEquals(delta[k].exp());
        }
    }

    /**
     * Perform the Inverse Laplace Transform.
     * @param function The function to apply the Inverse Laplace to.
     * @param time Argument at which to evaluate the time response.
     * @return {@code L<sup>-1</sup>{Y(s)} = y(t) evaluated at t = time.}
     */
    public double inverseTransform(ComplexUnivariateFunction function, double time) {
        if (time == 0.0 || time == -0.0) {
            time = ConstantsETK.DOUBLE_EPS;
        }
        double fb = 0.0;
        for(int k = 0; k < M; ++k) {
            fb += gamma[k].multiply(function.evaluateAt(delta[k].divide(time))).real();
        }
        return 0.4 / time * fb;
    }
}
