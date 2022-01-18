package com.wildbitsfoundry.etk4j.math.specialfunctions;

public final class Gamma {

    private Gamma() {
    }

    // TODO test

    /**
     * Gamma function
     * @param x Argument at which to evaluate the function.
     * @param odd Output: the odd part.
     * @param even Output: the even part.
     * @return the value of the gamma function evaluated at <code>x</code>.
     */
    public static double gamma(double x, double[] odd, double[] even) {
        boolean inv;
        double y, s, f = 0, g;

        if (x < 0.5) {
            y = x - Math.floor(x / 2.0) * 2;
            s = 3.14159265358979;
            if (y >= 1.0) {
                s = -s;
                y = 2.0 - y;
            }
            if (y >= 0.5) {
                y = 1.0 - y;
            }
            inv = true;
            x = 1.0 - x;
            f = s / Math.sin(3.14159265358979 * y);
        } else {
            inv = false;
        }
        if (x > 22.0)
            g = Math.exp(loggamma(x));
        else {
            s = 1.0;
            while (x > 1.5) {
                x = x - 1.0;
                s *= x;
            }
            g = s / recipgamma(1.0 - x, odd, even);
        }
        return (inv ? f / g : g);
    }

    // TODO test

    /**
     * Reciprocal of the gamma function.
     * @param x Argument at which to evaluate the function. x must be between [-0.5, 0.5].
     * @param odd Output: the odd part.
     * @param even Output: the even part.
     * @return the value of the reciprocal gamma function evaluated at 1 - x.
     */
    public static double recipgamma(double x, double[] odd, double[] even) {
        int i;
        double alfa, beta, x2;
        double[] b = new double[13];
        b[1] = -0.283876542276024;
        b[2] = -0.076852840844786;
        b[3] = 0.001706305071096;
        b[4] = 0.001271927136655;
        b[5] = 0.000076309597586;
        b[6] = -0.000004971736704;
        b[7] = -0.000000865920800;
        b[8] = -0.000000033126120;
        b[9] = 0.000000001745136;
        b[10] = 0.000000000242310;
        b[11] = 0.000000000009161;
        b[12] = -0.000000000000170;
        x2 = x * x * 8.0;
        alfa = -0.000000000000001;
        beta = 0.0;
        for (i = 12; i >= 2; i -= 2) {
            beta = -(alfa * 2.0 + beta);
            alfa = -beta * x2 - alfa + b[i];
        }
        even[0] = (beta / 2.0 + alfa) * x2 - alfa + 0.921870293650453;
        alfa = -0.000000000000034;
        beta = 0.0;
        for (i = 11; i >= 1; i -= 2) {
            beta = -(alfa * 2.0 + beta);
            alfa = -beta * x2 - alfa + b[i];
        }
        odd[0] = (alfa + beta) * 2.0;
        return odd[0] * x + even[0];
    }

    // TODO test
    /**
     * Natural logarithm of the gamma function.
     * @param x Argument at which to evaluate the function. x must be >= 0.
     * @return the natural logarithm of the gamma function at <code>x</code>.
     */
    public static double loggamma(double x) {
        int i;
        double r, x2, y, f, u0, u1, u, z;
        double[] b = new double[19];

        if (x > 13.0) {
            r = 1.0;
            while (x <= 22.0) {
                r /= x;
                x += 1.0;
            }
            x2 = -1.0 / (x * x);
            r = Math.log(r);
            return Math.log(x) * (x - 0.5) - x + r + 0.918938533204672 +
                    (((0.595238095238095e-3 * x2 + 0.793650793650794e-3) * x2 +
                            0.277777777777778e-2) * x2 + 0.833333333333333e-1) / x;
        } else {
            f = 1.0;
            u0 = u1 = 0.0;
            b[1] = -0.0761141616704358;
            b[2] = 0.0084323249659328;
            b[3] = -0.0010794937263286;
            b[4] = 0.0001490074800369;
            b[5] = -0.0000215123998886;
            b[6] = 0.0000031979329861;
            b[7] = -0.0000004851693012;
            b[8] = 0.0000000747148782;
            b[9] = -0.0000000116382967;
            b[10] = 0.0000000018294004;
            b[11] = -0.0000000002896918;
            b[12] = 0.0000000000461570;
            b[13] = -0.0000000000073928;
            b[14] = 0.0000000000011894;
            b[15] = -0.0000000000001921;
            b[16] = 0.0000000000000311;
            b[17] = -0.0000000000000051;
            b[18] = 0.0000000000000008;
            if (x < 1.0) {
                f = 1.0 / x;
                x += 1.0;
            } else
                while (x > 2.0) {
                    x -= 1.0;
                    f *= x;
                }
            f = Math.log(f);
            y = x + x - 3.0;
            z = y + y;
            for (i = 18; i >= 1; i--) {
                u = u0;
                u0 = z * u0 + b[i] - u1;
                u1 = u;
            }
            return (u0 * y + 0.491415393029387 - u1) * (x - 1.0) * (x - 2.0) + f;
        }
    }
}
