package com.wildbitsfoundry.etk4j.signals.laplace;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.signals.filters.ButterWorth;
import com.wildbitsfoundry.etk4j.signals.filters.Elliptic;
import com.wildbitsfoundry.etk4j.signals.filters.FilterOrderResults;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

import java.util.Arrays;

public class Laplace {

    public UnivariateFunction FunctionDelegate;
    static double[] V;       //  Stehfest coefficients
    static double ln2;       //  log of 2
    final static int DefaultStehfestN = 16;

    static {
        InitStehfest(DefaultStehfestN);
    }

    public static double Transform(UnivariateFunction F, double s) {
        final int DefaultIntegralN = 5000;
        double du = 0.5 / (double) DefaultIntegralN;
        double y = -F.evaluateAt(0) / 2.0;
        double u = 0;
        double limit = 1.0 - 1e-10;
        while (u < limit) {
            u += du;
            y += 2.0 * Math.pow(u, s - 1) * F.evaluateAt(-Math.log(u));
            u += du;
            y += Math.pow(u, s - 1) * F.evaluateAt(-Math.log(u));
        }
        return 2.0 * y * du / 3.0;
    }

    public static double InverseTransform(UnivariateFunction f, double t) {
        if (t == 0.0) {
            t = ConstantsETK.DOUBLE_EPS;
        } else if (t == -0.0) {
            t = -ConstantsETK.DOUBLE_EPS;
        }
        double ln2t = ln2 / t;
        double x = 0;
        double y = 0;
        for (int i = 0; i < V.length; i++) {
            x += ln2t;
            y += V[i] * f.evaluateAt(x);
        }
        return ln2t * y;
    }

    // http://www.columbia.edu/~ww2040/UnifiedDraft.pdf
    // TODO optimize this. All this coefficients can be calculated before hand
    public static double InverseTransformTalbot(RationalFunction f, double t, int M) {
        if (t == 0.0) {
            t = ConstantsETK.DOUBLE_EPS;
        } else if (t == -0.0) {
            t = -ConstantsETK.DOUBLE_EPS;
        }
        Complex[] delta = ComplexArrays.zeros(M);
        delta[0] = Complex.fromReal(2.0 * M / 5.0);
        for (int k = 1; k < M; ++k) {
            delta[k] = new Complex(1.0 / Math.tan(k * Math.PI / M), 1.0);
            delta[k].multiplyEquals(2.0 * k * Math.PI / 5.0);
        }
        Complex[] gamma = ComplexArrays.zeros(M);
        gamma[0] = delta[0].exp().multiply(0.5);
        for (int k = 1; k < M; ++k) {
            double cotTerm = 1.0 / Math.tan(k * Math.PI / M);
            gamma[k] = new Complex(1.0, (k * Math.PI / M) * (1.0 + cotTerm * cotTerm));
            gamma[k].subtractEquals(Complex.fromImaginary(cotTerm));
            gamma[k].multiplyEquals(delta[k].exp());
        }
        double fb = 0.0;
        for(int k = 0; k < M; ++k) {
            fb += gamma[k].multiply(f.evaluateAt(delta[k].divide(t))).real();
        }
        return 0.4 / t * fb;
    }

    public static double Factorial(int N) {
        double x = 1;
        if (N > 1) {
            for (int i = 2; i <= N; i++)
                x = i * x;
        }
        return x;
    }

    //public static double Integrate(FunctionDelegate f, double Min, double Max)
    //{
    //    return Integrate(f, Min, Max, 100);
    //}

    //public static double Integrate(FunctionDelegate f, double XMin, double XMax, int N)
    //{
    //    double dx = (XMax - XMin) / (double)N / 2.0;
    //    double y = (f(XMin) - f(XMax))/2.0;
    //    double x = XMin;
    //    double limit = XMax - 1e-10;
    //    while (x < limit)
    //    {
    //        x += dx;
    //        y += 2.0*f(x);
    //        x += dx;
    //        y += f(x);
    //    }
    //    return 2.0 * y * dx / 3.0;
    //}

    public static void InitStehfest() {
        InitStehfest(DefaultStehfestN);
    }

    public static void InitStehfest(int N) {
        ln2 = Math.log(2.0);
        int N2 = N / 2;
        int NV = 2 * N2;
        V = new double[NV];
        int sign = 1;
        if ((N2 % 2) != 0)
            sign = -1;
        for (int i = 0; i < NV; i++) {
            int kmin = (i + 2) / 2;
            int kmax = i + 1;
            if (kmax > N2)
                kmax = N2;
            V[i] = 0;
            sign = -sign;
            for (int k = kmin; k <= kmax; k++) {
                V[i] = V[i] + (Math.pow(k, N2) / Factorial(k)) * (Factorial(2 * k)
                        / Factorial(2 * k - i - 1)) / Factorial(N2 - k) / Factorial(k - 1)
                        / Factorial(i + 1 - k);
            }
            V[i] = sign * V[i];
        }
    }

    public static void main(String[] args) {
        RationalFunction rf = new RationalFunction(new double[]{1.0}, new double[]{1.0, 0.0});
        double[] t = {1, 10, 100, 1000, 10000};
        for (int i = 0; i < t.length; ++i) {
            t[i] = InverseTransformTalbot(rf, t[i], 64);
        }
        System.out.println(Arrays.toString(t));

        // Specs for low pass filter
        FilterSpecs.LowPassSpecs lpSpecs = new FilterSpecs.LowPassSpecs();
        lpSpecs.setPassBandRipple(1.5); // 1.5 dB gain/ripple refer to note
        lpSpecs.setStopBandAttenuation(60.0); // 60 dB at the stop band
        lpSpecs.setPassBandFrequency(2500); // 2500 Hz cutoff frequency
        lpSpecs.setStopBandFrequency(10000); // 10000 Hz stop band frequency

        FilterOrderResults.OrderAndCutoffFrequency nW0 = ButterWorth.buttord(lpSpecs);
        nW0 = Elliptic.ellipord(lpSpecs);
        TransferFunction el = Elliptic.newLowPass(nW0.getOrder(), lpSpecs.getPassBandRipple(),
                lpSpecs.getStopBandAttenuation(), nW0.getCutoffFrequency());

        double[] time = {0, 0.000124391246294398, 0.000248782492588796, 0.000373173738883193, 0.000497564985177591, 0.000621956231471989, 0.000746347477766387, 0.000870738724060784, 0.000995129970355182, 0.00111952121664958, 0.00124391246294398, 0.00136830370923838, 0.00149269495553277, 0.00161708620182717, 0.00174147744812157, 0.00186586869441597, 0.00199025994071036, 0.00211465118700476, 0.00223904243329916, 0.00236343367959356, 0.00248782492588796, 0.00261221617218235, 0.00273660741847675, 0.00286099866477115, 0.00298538991106555, 0.00310978115735994, 0.00323417240365434, 0.00335856364994874, 0.00348295489624314, 0.00360734614253753, 0.00373173738883193, 0.00385612863512633, 0.00398051988142073, 0.00410491112771513, 0.00422930237400952, 0.00435369362030392, 0.00447808486659832, 0.00460247611289272, 0.00472686735918712, 0.00485125860548151, 0.00497564985177591, 0.00510004109807031, 0.00522443234436471, 0.00534882359065910, 0.00547321483695350, 0.00559760608324790, 0.00572199732954230, 0.00584638857583669, 0.00597077982213109, 0.00609517106842549, 0.00621956231471989, 0.00634395356101429, 0.00646834480730868, 0.00659273605360308, 0.00671712729989748, 0.00684151854619188, 0.00696590979248627, 0.00709030103878067, 0.00721469228507507, 0.00733908353136947, 0.00746347477766387, 0.00758786602395826, 0.00771225727025266, 0.00783664851654706, 0.00796103976284146, 0.00808543100913585, 0.00820982225543025, 0.00833421350172465, 0.00845860474801905, 0.00858299599431345, 0.00870738724060784, 0.00883177848690224, 0.00895616973319664, 0.00908056097949104, 0.00920495222578543, 0.00932934347207983, 0.00945373471837423, 0.00957812596466863, 0.00970251721096303, 0.00982690845725742, 0.00995129970355182, 0.0100756909498462, 0.0102000821961406, 0.0103244734424350, 0.0104488646887294, 0.0105732559350238, 0.0106976471813182, 0.0108220384276126, 0.0109464296739070, 0.0110708209202014, 0.0111952121664958, 0.0113196034127902, 0.0114439946590846, 0.0115683859053790, 0.0116927771516734, 0.0118171683979678, 0.0119415596442622, 0.0120659508905566, 0.0121903421368510, 0.0123147333831454, 0.0124391246294398, 0.0125635158757342, 0.0126879071220286, 0.0128122983683230, 0.0129366896146174, 0.0130610808609118, 0.0131854721072062, 0.0133098633535006, 0.0134342545997950, 0.0135586458460894, 0.0136830370923838, 0.0138074283386782, 0.0139318195849725, 0.0140562108312669, 0.0141806020775613, 0.0143049933238557, 0.0144293845701501, 0.0145537758164445, 0.0146781670627389, 0.0148025583090333, 0.0149269495553277, 0.0150513408016221, 0.0151757320479165, 0.0153001232942109, 0.0154245145405053, 0.0155489057867997, 0.0156732970330941, 0.0157976882793885, 0.0159220795256829, 0.0160464707719773, 0.0161708620182717, 0.0162952532645661, 0.0164196445108605, 0.0165440357571549, 0.0166684270034493, 0.0167928182497437, 0.0169172094960381, 0.0170416007423325, 0.0171659919886269};
        double[] yout = el.step(time);
        System.out.println(Arrays.toString(yout));

    }
}
