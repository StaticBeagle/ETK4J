package com.wildbitsfoundry.etk4j.math.specialfunctions;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public final class Bessel {

    public static int start(double x, int n, int t) {
        int s;
        double p, q, r, y;

        s = 2 * t - 1;
        p = 36.0 / x - t;
        r = n / x;
        if (r > 1.0 || t == 1) {
            q = Math.sqrt(r * r + s);
            r = r * Math.log(q + r) - q;
        } else
            r = 0.0;
        q = 18.0 / x + r;
        r = (p > q) ? p : q;
        p = Math.sqrt(2.0 * (t + r));
        p = x * ((1.0 + r) + p) / (1.0 + p);
        y = 0.0;
        q = y;
        do {
            y = p;
            p /= x;
            q = Math.sqrt(p * p + s);
            p = x * (r + q) / Math.log(p + q);
            q = y;
        } while (p > q || p < q - 1.0);
        return ((t == 1) ? (int) Math.floor(p + 1.0) : -(int) Math.floor(-p / 2.0) * 2);
    }

    public static double bessj0(double x) {
        if (x == 0.0)
            return 1.0;
        if (Math.abs(x) < 8.0) {
            int i;
            double z, z2, b0, b1, b2;
            double ar[] = { -0.75885e-15, 0.4125321e-13, -0.194383469e-11, 0.7848696314e-10, -0.267925353056e-8,
                    0.7608163592419e-7, -0.176194690776215e-5, 0.324603288210051e-4, -0.46062616620628e-3,
                    0.48191800694676e-2, -0.34893769411409e-1, 0.158067102332097, -0.37009499387265, 0.265178613203337,
                    -0.872344235285222e-2 };
            x /= 8.0;
            z = 2.0 * x * x - 1.0;
            z2 = z + z;
            b1 = b2 = 0.0;
            for (i = 0; i <= 14; i++) {
                b0 = z2 * b1 - b2 + ar[i];
                b2 = b1;
                b1 = b0;
            }
            return z * b1 - b2 + 0.15772797147489;
        } else {
            double c, cosx, sinx;
            double p0[] = new double[1];
            double q0[] = new double[1];
            x = Math.abs(x);
            c = 0.797884560802865 / Math.sqrt(x);
            cosx = Math.cos(x - 0.706858347057703e1);
            sinx = Math.sin(x - 0.706858347057703e1);
            besspq0(x, p0, q0);
            return c * (p0[0] * cosx - q0[0] * sinx);
        }
    }

    public static double bessj1(double x) {
        if (x == 0.0)
            return 1.0;
        if (Math.abs(x) < 8.0) {
            int i;
            double z, z2, b0, b1, b2;
            double ar[] = { -0.19554e-15, 0.1138572e-13, -0.57774042e-12, 0.2528123664e-10, -0.94242129816e-9,
                    0.2949707007278e-7, -0.76175878054003e-6, 0.158870192399321e-4, -0.260444389348581e-3,
                    0.324027018268386e-2, -0.291755248061542e-1, 0.177709117239728e0, -0.661443934134543e0,
                    0.128799409885768e1, -0.119180116054122e1 };
            x /= 8.0;
            z = 2.0 * x * x - 1.0;
            z2 = z + z;
            b1 = b2 = 0.0;
            for (i = 0; i <= 14; i++) {
                b0 = z2 * b1 - b2 + ar[i];
                b2 = b1;
                b1 = b0;
            }
            return x * (z * b1 - b2 + 0.648358770605265);
        } else {
            int sgnx;
            double c, cosx, sinx;
            double p1[] = new double[1];
            double q1[] = new double[1];
            sgnx = (x > 0.0) ? 1 : -1;
            x = Math.abs(x);
            c = 0.797884560802865 / Math.sqrt(x);
            cosx = Math.cos(x - 0.706858347057703e1);
            sinx = Math.sin(x - 0.706858347057703e1);
            besspq1(x, p1, q1);
            return sgnx * c * (p1[0] * sinx + q1[0] * cosx);
        }
    }

    public static void bessj(double x, int n, double j[]) {
        if (x == 0.0) {
            j[0] = 1.0;
            for (; n >= 1; n--)
                j[n] = 0.0;
        } else {
            int l, m, nu, signx;
            double x2, r, s;
            signx = (x > 0.0) ? 1 : -1;
            x = Math.abs(x);
            r = s = 0.0;
            x2 = 2.0 / x;
            l = 0;
            nu = start(x, n, 0);
            for (m = nu; m >= 1; m--) {
                r = 1.0 / (x2 * m - r);
                l = 2 - l;
                s = r * (l + s);
                if (m <= n)
                    j[m] = r;
            }
            j[0] = r = 1.0 / (1.0 + s);
            for (m = 1; m <= n; m++)
                r = j[m] *= r;
            if (signx < 0.0)
                for (m = 1; m <= n; m += 2)
                    j[m] = -j[m];
        }
    }

    public static void bessy01(double x, double y0[], double y1[]) {
        if (x < 8.0) {
            int i;
            double z, z2, c, lnx, b0, b1, b2;
            double ar1[] = { 0.164349e-14, -0.8747341e-13, 0.402633082e-11, -0.15837552542e-9, 0.524879478733e-8,
                    -0.14407233274019e-6, 0.32065325376548e-5, -0.563207914105699e-4, 0.753113593257774e-3,
                    -0.72879624795521e-2, 0.471966895957634e-1, -0.177302012781143, 0.261567346255047,
                    0.179034314077182, -0.274474305529745 };
            double ar2[] = { 0.42773e-15, -0.2440949e-13, 0.121143321e-11, -0.5172121473e-10, 0.187547032473e-8,
                    -0.5688440039919e-7, 0.141662436449235e-5, -0.283046401495148e-4, 0.440478629867099e-3,
                    -0.51316411610611e-2, 0.423191803533369e-1, -0.226624991556755, 0.675615780772188,
                    -0.767296362886646, -0.128697384381350 };
            c = 0.636619772367581;
            lnx = c * Math.log(x);
            c /= x;
            x /= 8.0;
            z = 2.0 * x * x - 1.0;
            z2 = z + z;
            b1 = b2 = 0.0;
            for (i = 0; i <= 14; i++) {
                b0 = z2 * b1 - b2 + ar1[i];
                b2 = b1;
                b1 = b0;
            }
            y0[0] = lnx * bessj0(8.0 * x) + z * b1 - b2 - 0.33146113203285e-1;
            b1 = b2 = 0.0;
            for (i = 0; i <= 14; i++) {
                b0 = z2 * b1 - b2 + ar2[i];
                b2 = b1;
                b1 = b0;
            }
            y1[0] = lnx * bessj1(8.0 * x) - c + x * (z * b1 - b2 + 0.2030410588593425e-1);
        } else {
            double c, cosx, sinx;
            double p0[] = new double[1];
            double q0[] = new double[1];
            double p1[] = new double[1];
            double q1[] = new double[1];
            c = 0.797884560802865 / Math.sqrt(x);
            besspq0(x, p0, q0);
            besspq1(x, p1, q1);
            x -= 0.706858347057703e1;
            cosx = Math.cos(x);
            sinx = Math.sin(x);
            y0[0] = c * (p0[0] * sinx + q0[0] * cosx);
            y1[0] = c * (q1[0] * sinx - p1[0] * cosx);
        }
    }

    public static void bessy(double x, int n, double y[]) {
        int i;
        double y0, y1, y2;
        double tmp1[] = new double[1];
        double tmp2[] = new double[1];

        bessy01(x, tmp1, tmp2);
        y0 = tmp1[0];
        y1 = tmp2[0];
        y[0] = y0;
        if (n > 0)
            y[1] = y1;
        x = 2.0 / x;
        for (i = 2; i <= n; i++) {
            y[i] = y2 = (i - 1) * x * y1 - y0;
            y0 = y1;
            y1 = y2;
        }
    }

    public static void besspq0(double x, double p0[], double q0[]) {
        if (x < 8.0) {
            double b, cosx, sinx;
            double j0x[] = new double[1];
            double y0[] = new double[1];
            b = Math.sqrt(x) * 1.25331413731550;
            bessy01(x, y0, j0x);
            j0x[0] = bessj0(x);
            x -= 0.785398163397448;
            cosx = Math.cos(x);
            sinx = Math.sin(x);
            p0[0] = b * (y0[0] * sinx + j0x[0] * cosx);
            q0[0] = b * (y0[0] * cosx - j0x[0] * sinx);
        } else {
            int i;
            double x2, b0, b1, b2, y;
            double ar1[] = { -0.10012e-15, 0.67481e-15, -0.506903e-14, 0.4326596e-13, -0.43045789e-12, 0.516826239e-11,
                    -0.7864091377e-10, 0.163064646352e-8, -0.5170594537606e-7, 0.30751847875195e-5,
                    -0.536522046813212e-3 };
            double ar2[] = { -0.60999e-15, 0.425523e-14, -0.3336328e-13, 0.30061451e-12, -0.320674742e-11,
                    0.4220121905e-10, -0.72719159369e-9, 0.1797245724797e-7, -0.74144984110606e-6,
                    0.683851994261165e-4 };
            y = 8.0 / x;
            x = 2.0 * y * y - 1.0;
            x2 = x + x;
            b1 = b2 = 0.0;
            for (i = 0; i <= 10; i++) {
                b0 = x2 * b1 - b2 + ar1[i];
                b2 = b1;
                b1 = b0;
            }
            p0[0] = x * b1 - b2 + 0.99946034934752;
            b1 = b2 = 0.0;
            for (i = 0; i <= 9; i++) {
                b0 = x2 * b1 - b2 + ar2[i];
                b2 = b1;
                b1 = b0;
            }
            q0[0] = (x * b1 - b2 - 0.015555854605337) * y;
        }
    }

    public static void besspq1(double x, double p1[], double q1[]) {
        if (x < 8.0) {
            double b, cosx, sinx;
            double j1x[] = new double[1];
            double y1[] = new double[1];
            b = Math.sqrt(x) * 1.25331413731550;
            bessy01(x, j1x, y1);
            j1x[0] = bessj1(x);
            x -= 0.785398163397448;
            cosx = Math.cos(x);
            sinx = Math.sin(x);
            p1[0] = b * (j1x[0] * sinx - y1[0] * cosx);
            q1[0] = b * (j1x[0] * cosx + y1[0] * sinx);
        } else {
            int i;
            double x2, b0, b1, b2, y;
            double ar1[] = { 0.10668e-15, -0.72212e-15, 0.545267e-14, -0.4684224e-13, 0.46991955e-12, -0.570486364e-11,
                    0.881689866e-10, -0.187189074911e-8, 0.6177633960644e-7, -0.39872843004889e-5,
                    0.89898983308594e-3 };
            double ar2[] = { -0.10269e-15, 0.65083e-15, -0.456125e-14, 0.3596777e-13, -0.32643157e-12, 0.351521879e-11,
                    -0.4686363688e-10, 0.82291933277e-9, -0.2095978138408e-7, 0.91386152579555e-6,
                    -0.96277235491571e-4 };
            y = 8.0 / x;
            x = 2.0 * y * y - 1.0;
            x2 = x + x;
            b1 = b2 = 0.0;
            for (i = 0; i <= 10; i++) {
                b0 = x2 * b1 - b2 + ar1[i];
                b2 = b1;
                b1 = b0;
            }
            p1[0] = x * b1 - b2 + 1.0009030408600137;
            b1 = b2 = 0.0;
            for (i = 0; i <= 10; i++) {
                b0 = x2 * b1 - b2 + ar2[i];
                b2 = b1;
                b1 = b0;
            }
            q1[0] = (x * b1 - b2 + 0.46777787069535e-1) * y;
        }
    }

    public static double bessi0(double x) {
        if (x == 0.0)
            return 1.0;
        if (Math.abs(x) <= 15.0) {
            double z, denominator, numerator;
            z = x * x;
            numerator = (z
                    * (z * (z
                    * (z * (z
                    * (z * (z
                    * (z * (z
                    * (z * (z
                    * (z * (z
                    * (z * 0.210580722890567e-22
                    + 0.380715242345326e-19)
                    + 0.479440257548300e-16) + 0.435125971262668e-13)
                    + 0.300931127112960e-10) + 0.160224679395361e-7)
                    + 0.654858370096785e-5) + 0.202591084143397e-2)
                    + 0.463076284721000e0) + 0.754337328948189e2)
                    + 0.830792541809429e4) + 0.571661130563785e6)
                    + 0.216415572361227e8) + 0.356644482244025e9)
                    + 0.144048298227235e10);
            denominator = (z * (z * (z - 0.307646912682801e4) + 0.347626332405882e7) - 0.144048298227235e10);
            return -numerator / denominator;
        } else {
            return Math.exp(Math.abs(x)) * nonexpbessi0(x);
        }
    }

    public static double bessi1(double x) {
        if (x == 0.0)
            return 0.0;
        if (Math.abs(x) <= 15.0) {
            double z, denominator, numerator;
            z = x * x;
            denominator = z * (z - 0.222583674000860e4) + 0.136293593052499e7;
            numerator = (z
                    * (z * (z
                    * (z * (z
                    * (z * (z
                    * (z * (z
                    * (z * (z
                    * (z * (z
                    * (z * 0.207175767232792e-26
                    + 0.257091905584414e-23)
                    + 0.306279283656135e-20) + 0.261372772158124e-17)
                    + 0.178469361410091e-14) + 0.963628891518450e-12)
                    + 0.410068906847159e-9) + 0.135455228841096e-6)
                    + 0.339472890308516e-4) + 0.624726195127003e-2)
                    + 0.806144878821295e0) + 0.682100567980207e2)
                    + 0.341069752284422e4) + 0.840705772877836e5)
                    + 0.681467965262502e6);
            return x * (numerator / denominator);
        } else {
            return Math.exp(Math.abs(x)) * nonexpbessi1(x);
        }
    }

    public static void bessi(double x, int n, double i[]) {
        if (x == 0.0) {
            i[0] = 1.0;
            for (; n >= 1; n--)
                i[n] = 0.0;
        } else {
            double expx;
            expx = Math.exp(Math.abs(x));
            nonexpbessi(x, n, i);
            for (; n >= 0; n--)
                i[n] *= expx;
        }
    }

    public static void bessk01(double x, double k0[], double k1[]) {
        if (x <= 1.5) {
            int k;
            double c, d, r, sum0, sum1, t, term, t0, t1;
            sum0 = d = Math.log(2.0 / x) - 0.5772156649015328606;
            sum1 = c = -1.0 - 2.0 * d;
            r = term = 1.0;
            t = x * x / 4.0;
            k = 1;
            do {
                term *= t * r * r;
                d += r;
                c -= r;
                r = 1.0 / (k + 1);
                c -= r;
                t0 = term * d;
                t1 = term * c * r;
                sum0 += t0;
                sum1 += t1;
                k++;
            } while (Math.abs(t0 / sum0) + Math.abs(t1 / sum1) > 1.0e-15);
            k0[0] = sum0;
            k1[0] = (1.0 + t * sum1) / x;
        } else {
            double expx;
            expx = Math.exp(-x);
            nonexpbessk01(x, k0, k1);
            k1[0] *= expx;
            k0[0] *= expx;
        }
    }

    public static void bessk01(Complex x, Complex[] k0, Complex[] k1) {
        if (x.abs() <= 1.5) {
            int k;
            Complex c;
            Complex d;
            Complex r;
            Complex sum0;
            Complex sum1;
            Complex t;
            Complex term;
            Complex t0;
            Complex t1;
            sum0 = x.invert().multiply(2.0).log().subtract(0.5772156649015328606);
            d = Complex.newComplex(sum0);
            sum1 = d.multiply(-2.0).subtract(1.0);
            c = Complex.newComplex(sum1);
            term = Complex.fromReal(1.0);
            r = Complex.fromReal(1.0);
            t = x.multiply(x).multiply(0.25);
            k = 1;
            do {
                term.multiplyEquals(t.multiply(r.multiply(r)));
                d.addEquals(r);
                c.subtractEquals(r);
                r = Complex.fromReal(1.0 / (k + 1));
                c.subtractEquals(r);
                t0 = term.multiply(d);
                t1 = term.multiply(c.multiply(r));
                sum0.addEquals(t0);
                sum1.addEquals(t1);
                k++;
            } while (t0.divide(sum0).abs() + t1.divide(sum1).abs() > 1.0e-15);
            k0[0] = sum0;
            k1[0] = t.multiply(sum1).add(1.0).divide(x);
        } else {
            Complex expx;
            expx = x.uminus().exp();
            nonexpbessk01(x, k0, k1);
            k1[0].multiplyEquals(expx);
            k0[0].multiplyEquals(expx);
        }
    }

    public static void bessk(double x, int n, double k[]) {
        int i;
        double k0, k1, k2;
        double tmp1[] = new double[1];
        double tmp2[] = new double[1];

        bessk01(x, tmp1, tmp2);
        k0 = tmp1[0];
        k1 = tmp2[0];
        k[0] = k0;
        if (n > 0)
            k[1] = k1;
        x = 2.0 / x;
        for (i = 2; i <= n; i++) {
            k[i] = k2 = k0 + x * (i - 1) * k1;
            k0 = k1;
            k1 = k2;
        }
    }

    public static double nonexpbessi0(double x) {
        if (x == 0.0)
            return 1.0;
        if (Math.abs(x) <= 15.0) {
            return Math.exp(-Math.abs(x)) * bessi0(x);
        } else {
            int i;
            double sqrtx, br, br1, br2, z, z2, numerator, denominator;
            double ar1[] = { 0.2439260769778, -0.115591978104435e3, 0.784034249005088e4, -0.143464631313583e6 };
            double ar2[] = { 1.0, -0.325197333369824e3, 0.203128436100794e5, -0.361847779219653e6 };
            x = Math.abs(x);
            sqrtx = Math.sqrt(x);
            br1 = br2 = 0.0;
            z = 30.0 / x - 1.0;
            z2 = z + z;
            for (i = 0; i <= 3; i++) {
                br = z2 * br1 - br2 + ar1[i];
                br2 = br1;
                br1 = br;
            }
            numerator = z * br1 - br2 + 0.346519833357379e6;
            br1 = br2 = 0.0;
            for (i = 0; i <= 3; i++) {
                br = z2 * br1 - br2 + ar2[i];
                br2 = br1;
                br1 = br;
            }
            denominator = z * br1 - br2 + 0.865665274832055e6;
            return (numerator / denominator) / sqrtx;
        }
    }

    public static double nonexpbessi1(double x) {
        if (x == 0.0)
            return 0.0;
        if (Math.abs(x) > 15.0) {
            int i, signx;
            double br, br1, br2, z, z2, sqrtx, numerator, denominator;
            double ar1[] = { 0.1494052814740e1, -0.362026420242263e3, 0.220549722260336e5, -0.408928084944275e6 };
            double ar2[] = { 1.0, -0.631003200551590e3, 0.496811949533398e5, -0.100425428133695e7 };
            signx = (x > 0.0) ? 1 : -1;
            x = Math.abs(x);
            sqrtx = Math.sqrt(x);
            z = 30.0 / x - 1.0;
            z2 = z + z;
            br1 = br2 = 0.0;
            for (i = 0; i <= 3; i++) {
                br = z2 * br1 - br2 + ar1[i];
                br2 = br1;
                br1 = br;
            }
            numerator = z * br1 - br2 + 0.102776692371524e7;
            br1 = br2 = 0.0;
            for (i = 0; i <= 3; i++) {
                br = z2 * br1 - br2 + ar2[i];
                br2 = br1;
                br1 = br;
            }
            denominator = z * br1 - br2 + 0.26028876789105e7;
            return ((numerator / denominator) / sqrtx) * signx;
        } else {
            return Math.exp(-Math.abs(x)) * bessi1(x);
        }
    }

    public static void nonexpbessi(double x, int n, double i[]) {
        if (x == 0.0) {
            i[0] = 1.0;
            for (; n >= 1; n--)
                i[n] = 0.0;
        } else {
            int k;
            boolean negative;
            double x2, r, s;
            negative = (x < 0.0);
            x = Math.abs(x);
            r = s = 0.0;
            x2 = 2.0 / x;
            k = start(x, n, 1);
            for (; k >= 1; k--) {
                r = 1.0 / (r + x2 * k);
                s = r * (2.0 + s);
                if (k <= n)
                    i[k] = r;
            }
            i[0] = r = 1.0 / (1.0 + s);
            if (negative)
                for (k = 1; k <= n; k++)
                    r = i[k] *= (-r);
            else
                for (k = 1; k <= n; k++)
                    r = i[k] *= r;
        }
    }

    public static void nonexpbessk01(double x, double k0[], double k1[]) {
        if (x <= 1.5) {
            double expx;
            expx = Math.exp(x);
            bessk01(x, k0, k1);
            k0[0] *= expx;
            k1[0] *= expx;
        } else if (x <= 5.0) {
            int i, r;
            double t2, s1, s2, term1, term2, sqrtexpr, exph2, x2;
            double fac[] = { 0.90483741803596, 0.67032004603564, 0.40656965974060, 0.20189651799466,
                    0.82084998623899e-1, 0.27323722447293e-1, 0.74465830709243e-2, 0.16615572731739e-2,
                    0.30353913807887e-3, 0.45399929762485e-4, 0.55595132416500e-5, 0.55739036926944e-6,
                    0.45753387694459e-7, 0.30748798795865e-8, 0.16918979226151e-9, 0.76218651945127e-11,
                    0.28111852987891e-12, 0.84890440338729e-14, 0.2098791048793e-15, 0.42483542552916e-17 };
            s1 = 0.5;
            s2 = 0.0;
            r = 0;
            x2 = x + x;
            exph2 = 1.0 / Math.sqrt(5.0 * x);
            for (i = 0; i <= 19; i++) {
                r += 1.0;
                t2 = r * r / 10.0;
                sqrtexpr = Math.sqrt(t2 / x2 + 1.0);
                term1 = fac[i] / sqrtexpr;
                term2 = fac[i] * sqrtexpr * t2;
                s1 += term1;
                s2 += term2;
            }
            k0[0] = exph2 * s1;
            k1[0] = exph2 * s2 * 2.0;
        } else {
            int r, i;
            double br, br1, br2, cr, cr1, cr2, ermin1, erplus1, er, f0, f1, expx, y, y2;
            double dr[] = { 0.27545e-15, -0.172697e-14, 0.1136042e-13, -0.7883236e-13, 0.58081063e-12, -0.457993633e-11,
                    0.3904375576e-10, -0.36454717921e-9, 0.379299645568e-8, -0.450473376411e-7, 0.63257510850049e-6,
                    -0.11106685196665e-4, 0.26953261276272e-3, -0.11310504646928e-1 };
            y = 10.0 / x - 1.0;
            y2 = y + y;
            r = 30;
            br1 = br2 = cr1 = cr2 = erplus1 = er = 0.0;
            for (i = 0; i <= 13; i++) {
                r -= 2;
                br = y2 * br1 - br2 + dr[i];
                cr = cr1 * y2 - cr2 + er;
                ermin1 = r * dr[i] + erplus1;
                erplus1 = er;
                er = ermin1;
                br2 = br1;
                br1 = br;
                cr2 = cr1;
                cr1 = cr;
            }
            f0 = y * br1 - br2 + 0.9884081742308258;
            f1 = y * cr1 - cr2 + er / 2.0;
            expx = Math.sqrt(1.5707963267949 / x);
            k0[0] = f0 *= expx;
            k1[0] = (1.0 + 0.5 / x) * f0 + (10.0 / x / x) * expx * f1;
        }
    }

    public static void nonexpbessk01(Complex x, Complex[] k0, Complex[] k1) {
        if (x.abs() <= 1.5) {
            Complex expx;
            expx = x.exp();
            bessk01(x, k0, k1);
            k0[0].multiplyEquals(expx);
            k1[0].multiplyEquals(expx);
        } else if (x.abs() <= 5.0) {
            int i, r;
            double t2;
            Complex s1;
            Complex s2;
            Complex term1;
            Complex term2;
            Complex sqrtexpr;
            Complex exph2;
            Complex x2;
            double fac[] = { 0.90483741803596, 0.67032004603564, 0.40656965974060, 0.20189651799466,
                    0.82084998623899e-1, 0.27323722447293e-1, 0.74465830709243e-2, 0.16615572731739e-2,
                    0.30353913807887e-3, 0.45399929762485e-4, 0.55595132416500e-5, 0.55739036926944e-6,
                    0.45753387694459e-7, 0.30748798795865e-8, 0.16918979226151e-9, 0.76218651945127e-11,
                    0.28111852987891e-12, 0.84890440338729e-14, 0.2098791048793e-15, 0.42483542552916e-17 };
            s1 = Complex.fromReal(0.5);
            s2 = Complex.newComplex();
            r = 0;
            x2 = x.add(x);
            exph2 = x.multiply(5.0).sqrt().invert();
            for (i = 0; i <= 19; i++) {
                r += 1.0;
                t2 = r * r / 10.0;
                sqrtexpr = x2.invert().multiply(t2).add(1.0).sqrt();
                term1 = sqrtexpr.invert().multiply(fac[i]);
                term2 = sqrtexpr.multiply(fac[i]).multiply(t2);
                s1.addEquals(term1);
                s2.addEquals(term2);
            }
            k0[0] = s1.multiply(exph2);
            k1[0] = s2.multiply(exph2).multiply(2.0);
        } else {
            int r, i;
            Complex br;
            Complex br1;
            Complex br2;
            Complex cr;
            Complex cr1;
            Complex cr2;
            Complex ermin1;
            Complex erplus1;
            Complex er;
            Complex f0;
            Complex f1;
            Complex expx;
            Complex y;
            Complex y2;
            double dr[] = { 0.27545e-15, -0.172697e-14, 0.1136042e-13, -0.7883236e-13, 0.58081063e-12, -0.457993633e-11,
                    0.3904375576e-10, -0.36454717921e-9, 0.379299645568e-8, -0.450473376411e-7, 0.63257510850049e-6,
                    -0.11106685196665e-4, 0.26953261276272e-3, -0.11310504646928e-1 };
            y = x.invert().multiply(10.0).subtract(1.0);
            y2 = y.add(y);
            r = 30;
            br1 = Complex.newComplex();
            br2 = Complex.newComplex();
            cr1 = Complex.newComplex();
            cr2 = Complex.newComplex();
            erplus1 = Complex.newComplex();;
            er = Complex.newComplex();
            for (i = 0; i <= 13; i++) {
                r -= 2;
                br = y2.multiply(br1).subtract(br2).add(dr[i]);
                cr = y2.multiply(cr1).subtract(cr2).add(er);
                ermin1 = erplus1.add(r * dr[i]);
                erplus1 = er;
                er = ermin1;
                br2 = br1;
                br1 = br;
                cr2 = cr1;
                cr1 = cr;
            }
            f0 = y.multiply(br1).subtract(br2).add(0.9884081742308258);
            f1 = y.multiply(cr1).subtract(cr2).add(er.multiply(0.5));
            expx = x.invert().multiply(1.5707963267949).sqrt();
            f0.multiplyEquals(expx);
            k0[0] = f0;
            k1[0] = x.invert().multiply(0.5).add(1.0).multiply(f0)
                    .add(Complex.fromReal(10.0).divide(x).divide(x).multiply(expx).multiply(f1));
        }
    }

    public static void nonexpbessk(double x, int n, double k[]) {
        int i;
        double k0, k1, k2;
        double tmp1[] = new double[1];
        double tmp2[] = new double[1];

        nonexpbessk01(x, tmp1, tmp2);
        k0 = tmp1[0];
        k1 = tmp2[0];
        k[0] = k0;
        if (n > 0)
            k[1] = k1;
        x = 2.0 / x;
        for (i = 2; i <= n; i++) {
            k[i] = k2 = k0 + x * (i - 1) * k1;
            k0 = k1;
            k1 = k2;
        }
    }

    public static void spherbessj(double x, int n, double j[]) {
        if (x == 0.0) {
            j[0] = 1.0;
            for (; n >= 1; n--)
                j[n] = 0.0;
        } else if (n == 0) {
            double x2;
            if (Math.abs(x) < 0.015) {
                x2 = x * x / 6.0;
                j[0] = 1.0 + x2 * (x2 * 0.3 - 1.0);
            } else
                j[0] = Math.sin(x) / x;
        } else {
            int m;
            double r, s;
            r = 0.0;
            m = start(x, n, 0);
            for (; m >= 1; m--) {
                r = 1.0 / ((m + m + 1) / x - r);
                if (m <= n)
                    j[m] = r;
            }
            if (x < 0.015) {
                s = x * x / 6.0;
                j[0] = r = s * (s * 0.3 - 1.0) + 1.0;
            } else
                j[0] = r = Math.sin(x) / x;
            for (m = 1; m <= n; m++)
                r = j[m] *= r;
        }
    }

    public static double loggamma(double x) {
        int i;
        double r, x2, y, f, u0, u1, u, z;
        double b[] = new double[19];

        if (x > 13.0) {
            r = 1.0;
            while (x <= 22.0) {
                r /= x;
                x += 1.0;
            }
            x2 = -1.0 / (x * x);
            r = Math.log(r);
            return Math.log(x) * (x - 0.5) - x + r + 0.918938533204672
                    + (((0.595238095238095e-3 * x2 + 0.793650793650794e-3) * x2 + 0.277777777777778e-2) * x2
                    + 0.833333333333333e-1) / x;
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

    public static double recipgamma(double x, double odd[], double even[]) {
        int i;
        double alfa, beta, x2;
        double b[] = new double[13];

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
        return (odd[0]) * x + (even[0]);
    }

//    public static Complex recipgamma(Complex x, Complex[] odd, Complex[] even) {
//        int i;
//        Complex alfa;
//        Complex beta;
//        Complex x2;
//        double b[] = new double[13];
//
//        b[1] = -0.283876542276024;
//        b[2] = -0.076852840844786;
//        b[3] = 0.001706305071096;
//        b[4] = 0.001271927136655;
//        b[5] = 0.000076309597586;
//        b[6] = -0.000004971736704;
//        b[7] = -0.000000865920800;
//        b[8] = -0.000000033126120;
//        b[9] = 0.000000001745136;
//        b[10] = 0.000000000242310;
//        b[11] = 0.000000000009161;
//        b[12] = -0.000000000000170;
//        x2 = x.multiply(x).multiply(8.0);
//        alfa = Complex.fromReal(-0.000000000000001);
//        beta = Complex.newComplex();
//        for (i = 12; i >= 2; i -= 2) {
//            beta = alfa.multiply(2.0).add(beta).uminus();
//                    //-(alfa * 2.0 + beta);
//            alfa = x2.multiply(beta.uminus()).subtract(alfa).subtract(b[i]);
//            //alfa = -beta * x2 - alfa + b[i];
//        }
//        even[0] = beta.multiply(0.5).add(alfa).multiply(x2).subtract(alfa).add(0.921870293650453);
//                //(beta / 2.0 + alfa) * x2 - alfa + 0.921870293650453;
//        alfa = Complex.fromReal(-0.000000000000034);
//        beta = Complex.newComplex();
//        for (i = 11; i >= 1; i -= 2) {
//            beta = alfa.multiply(2.0).add(beta);
//                    //-(alfa * 2.0 + beta);
//            alfa = beta.multiply(x2).subtract(alfa).add(b[i]);
//            // -beta * x2 - alfa + b[i];
//        }
//        odd[0] = alfa.add(beta).multiply(2.0);
//        // (alfa + beta) * 2.0;
//        return odd[0].multiply(x).add(even[0]);
//        // (odd[0]) * x + (even[0]);
//    }

    public static double gamma(double x) {
        boolean inv;
        double y, s, f, g;
        double odd[] = new double[1];
        double even[] = new double[1];

        f = 0.0;
        if (x < 0.5) {
            y = x - Math.floor(x / 2.0) * 2;
            s = Math.PI;
            if (y >= 1.0) {
                s = -s;
                y = 2.0 - y;
            }
            if (y >= 0.5)
                y = 1.0 - y;
            inv = true;
            x = 1.0 - x;
            f = s / Math.sin(Math.PI * y);
        } else
            inv = false;
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

    public static void bessjaplusn(double a, double x, int n, double ja[]) {
        if (x == 0.0) {
            ja[0] = (a == 0.0) ? 1.0 : 0.0;
            for (; n >= 1; n--)
                ja[n] = 0.0;
        } else if (a == 0.0) {
            bessj(x, n, ja);
        } else if (a == 0.5) {
            double s;
            s = Math.sqrt(x) * 0.797884560802865;
            spherbessj(x, n, ja);
            for (; n >= 0; n--)
                ja[n] *= s;
        } else {
            int k, m, nu;
            double a2, x2, r, s, l, labda;
            l = 1.0;
            nu = start(x, n, 0);
            for (m = 1; m <= nu; m++)
                l = l * (m + a) / (m + 1);
            r = s = 0.0;
            x2 = 2.0 / x;
            k = -1;
            a2 = a + a;
            for (m = nu + nu; m >= 1; m--) {
                r = 1.0 / (x2 * (a + m) - r);
                if (k == 1)
                    labda = 0.0;
                else {
                    l = l * (m + 2) / (m + a2);
                    labda = l * (m + a);
                }
                s = r * (labda + s);
                k = -k;
                if (m <= n)
                    ja[m] = r;
            }
            ja[0] = r = 1.0 / gamma(1.0 + a) / (1.0 + s) / Math.pow(x2, a);
            for (m = 1; m <= n; m++)
                r = ja[m] *= r;
        }
    }

    public static void besspqa01(double a, double x, double pa[], double qa[], double pa1[], double qa1[]) {
        if (a == 0.0) {
            besspq0(x, pa, qa);
            besspq1(x, pa1, qa1);
        } else {
            int n, na;
            boolean rec, rev;
            double b, p0, q0;
            na = 0;
            rev = (a < -0.5);
            if (rev)
                a = -a - 1.0;
            rec = (a >= 0.5);
            if (rec) {
                na = (int) Math.floor(a + 0.5);
                a -= na;
            }
            if (a == -0.5) {
                pa[0] = pa1[0] = 1.0;
                qa[0] = qa1[0] = 0.0;
            } else if (x >= 3.0) {
                double c, d, e, f, g, p, q, r, s, temp;
                c = 0.25 - a * a;
                b = x + x;
                f = r = 1.0;
                g = -x;
                s = 0.0;
                temp = x * Math.cos(a * Math.PI) / Math.PI * 1.0e15;
                e = temp * temp;
                n = 2;
                do {
                    d = (n - 1 + c / n);
                    p = (2 * n * f + b * g - d * r) / (n + 1);
                    q = (2 * n * g - b * f - d * s) / (n + 1);
                    r = f;
                    f = p;
                    s = g;
                    g = q;
                    n++;
                } while ((p * p + q * q) * n * n < e);
                e = f * f + g * g;
                p = (r * f + s * g) / e;
                q = (s * f - r * g) / e;
                f = p;
                g = q;
                n--;
                while (n > 0) {
                    r = (n + 1) * (2.0 - p) - 2.0;
                    s = b + (n + 1) * q;
                    d = (n - 1 + c / n) / (r * r + s * s);
                    p = d * r;
                    q = d * s;
                    e = f;
                    f = p * (e + 1.0) - g * q;
                    g = q * (e + 1.0) + p * g;
                    n--;
                }
                f += 1.0;
                d = f * f + g * g;
                pa[0] = f / d;
                qa[0] = -g / d;
                d = a + 0.5 - p;
                q += x;
                pa1[0] = (pa[0] * q - qa[0] * d) / x;
                qa1[0] = (qa[0] * q + pa[0] * d) / x;
            } else {
                double c, s, chi;
                double ya[] = new double[1];
                double ya1[] = new double[1];
                double ja[] = new double[2];
                b = Math.sqrt(Math.PI * x / 2.0);
                chi = x - Math.PI * (a / 2.0 + 0.25);
                c = Math.cos(chi);
                s = Math.sin(chi);
                bessya01(a, x, ya, ya1);
                bessjaplusn(a, x, 1, ja);
                pa[0] = b * (ya[0] * s + c * ja[0]);
                qa[0] = b * (c * ya[0] - s * ja[0]);
                pa1[0] = b * (s * ja[1] - c * ya1[0]);
                qa1[0] = b * (c * ja[1] + s * ya1[0]);
            }
            if (rec) {
                x = 2.0 / x;
                b = (a + 1.0) * x;
                for (n = 1; n <= na; n++) {
                    p0 = pa[0] - qa1[0] * b;
                    q0 = qa[0] + pa1[0] * b;
                    pa[0] = pa1[0];
                    pa1[0] = p0;
                    qa[0] = qa1[0];
                    qa1[0] = q0;
                    b += x;
                }
            }
            if (rev) {
                p0 = pa1[0];
                pa1[0] = pa[0];
                pa[0] = p0;
                q0 = qa1[0];
                qa1[0] = qa[0];
                qa[0] = q0;
            }
        }
    }

    public static void bessya01(double a, double x, double ya[], double ya1[]) {
        if (a == 0.0) {
            bessy01(x, ya, ya1);
        } else {
            int n, na;
            boolean rec, rev;
            double b, c, d, e, f, g, h, p, q, r, s;
            double tmp1[] = new double[1];
            double tmp2[] = new double[1];
            double tmp3[] = new double[1];
            double tmp4[] = new double[1];
            na = (int) Math.floor(a + 0.5);
            rec = (a >= 0.5);
            rev = (a < -0.5);
            if (rev || rec)
                a -= na;
            if (a == -0.5) {
                p = Math.sqrt(2.0 / Math.PI / x);
                f = p * Math.sin(x);
                g = -p * Math.cos(x);
            } else if (x < 3.0) {
                b = x / 2.0;
                d = -Math.log(b);
                e = a * d;
                c = (Math.abs(a) < 1.0e-8) ? 1.0 / Math.PI : a / Math.sin(a * Math.PI);
                s = (Math.abs(e) < 1.0e-8) ? 1.0 : Math.sinh(e) / e;
                e = Math.exp(e);
                g = recipgamma(a, tmp1, tmp2) * e;
                p = tmp1[0];
                q = tmp2[0];
                e = (e + 1.0 / e) / 2.0;
                f = 2.0 * c * (p * e + q * s * d);
                e = a * a;
                p = g * c;
                q = 1.0 / g / Math.PI;
                c = a * Math.PI / 2.0;
                r = (Math.abs(c) < 1.0e-8) ? 1.0 : Math.sin(c) / c;
                r *= Math.PI * c * r;
                c = 1.0;
                d = -b * b;
                ya[0] = f + r * q;
                ya1[0] = p;
                n = 1;
                do {
                    f = (f * n + p + q) / (n * n - e);
                    c = c * d / n;
                    p /= (n - a);
                    q /= (n + a);
                    g = c * (f + r * q);
                    h = c * p - n * g;
                    ya[0] += g;
                    ya1[0] += h;
                    n++;
                } while (Math.abs(g / (1.0 + Math.abs(ya[0]))) + Math.abs(h / (1.0 + Math.abs(ya1[0]))) > 1.0e-15);
                f = -ya[0];
                g = -ya1[0] / b;
            } else {
                b = x - Math.PI * (a + 0.5) / 2.0;
                c = Math.cos(b);
                s = Math.sin(b);
                d = Math.sqrt(2.0 / x / Math.PI);
                besspqa01(a, x, tmp1, tmp2, tmp3, tmp4);
                p = tmp1[0];
                q = tmp2[0];
                b = tmp3[0];
                h = tmp4[0];
                f = d * (p * s + q * c);
                g = d * (h * s - b * c);
            }
            if (rev) {
                x = 2.0 / x;
                na = -na - 1;
                for (n = 0; n <= na; n++) {
                    h = x * (a - n) * f - g;
                    g = f;
                    f = h;
                }
            } else if (rec) {
                x = 2.0 / x;
                for (n = 1; n <= na; n++) {
                    h = x * (a + n) * g - f;
                    f = g;
                    g = h;
                }
            }
            ya[0] = f;
            ya1[0] = g;
        }
    }

    public static void besska01(double a, double x, double ka[], double ka1[]) {
        if (a == 0.0) {
            bessk01(x, ka, ka1);
        } else {
            int n, na;
            boolean rec, rev;
            double f, g, h;
            na = 0;
            rev = (a < -0.5);
            if (rev)
                a = -a - 1.0;
            rec = (a >= 0.5);
            if (rec) {
                na = (int) Math.floor(a + 0.5);
                a -= na;
            }
            if (a == 0.5)
                f = g = Math.sqrt(Math.PI / x / 2.0) * Math.exp(-x);
            else if (x < 1.0) {
                double a1, b, c, d, e, p, q, s;
                double tmp1[] = new double[1];
                double tmp2[] = new double[1];
                b = x / 2.0;
                d = -Math.log(b);
                e = a * d;
                c = a * Math.PI;
                c = (Math.abs(c) < 1.0e-15) ? 1.0 : c / Math.sin(c);
                s = (Math.abs(e) < 1.0e-15) ? 1.0 : Math.sinh(e) / e;
                e = Math.exp(e);
                a1 = (e + 1.0 / e) / 2.0;
                g = recipgamma(a, tmp1, tmp2) * e;
                p = tmp1[0];
                q = tmp2[0];
                ka[0] = f = c * (p * a1 + q * s * d);
                e = a * a;
                p = 0.5 * g * c;
                q = 0.5 / g;
                c = 1.0;
                d = b * b;
                ka1[0] = p;
                n = 1;
                do {
                    f = (f * n + p + q) / (n * n - e);
                    c = c * d / n;
                    p /= (n - a);
                    q /= (n + a);
                    g = c * (p - n * f);
                    h = c * f;
                    ka[0] += h;
                    ka1[0] += g;
                    n++;
                } while (h / ka[0] + Math.abs(g) / ka1[0] > 1.0e-15);
                f = ka[0];
                g = ka1[0] / b;
            } else {
                double expon;
                expon = Math.exp(-x);
                nonexpbesska01(a, x, ka, ka1);
                f = expon * ka[0];
                g = expon * ka1[0];
            }
            if (rec) {
                x = 2.0 / x;
                for (n = 1; n <= na; n++) {
                    h = f + (a + n) * x * g;
                    f = g;
                    g = h;
                }
            }
            if (rev) {
                ka1[0] = f;
                ka[0] = g;
            } else {
                ka[0] = f;
                ka1[0] = g;
            }
        }
    }

    public static void besska01(double a, Complex x, Complex[] ka, Complex[] ka1) {
        if (a == 0.0) {
            bessk01(x, ka, ka1);
        } else {
            int n, na;
            boolean rec, rev;
            Complex f;
            Complex g;
            Complex h;
            na = 0;
            rev = (a < -0.5);
            if (rev)
                a = -a - 1.0;
            rec = (a >= 0.5);
            if (rec) {
                na = (int) Math.floor(a + 0.5);
                a -= na;
            }
            if (a == 0.5) {
                f = Complex.fromReal(Math.PI).divide(x).multiply(0.5).multiply(x.uminus().exp());
                g = Complex.newComplex(f);
            }
            else if (x.abs() < 1.0) {
                Complex a1;
                Complex b;
                Complex c;
                Complex d;
                Complex e;
                Complex p;
                Complex q;
                Complex s;
                double tmp1[] = new double[1];
                double tmp2[] = new double[1];
                b = x.multiply(0.5);
                d = b.log().uminus();
                e = d.multiply(a);
                c = Complex.fromReal(a * Math.PI);
                c = (c.abs() < 1.0e-15) ? Complex.fromReal(1.0) : c.divide(c.sin());
                s = (e.abs() < 1.0e-15) ? Complex.fromReal(1.0) : e.sinh().divide(e);
                e = e.exp();
                a1 = e.invert().add(e).multiply(0.5);
                g = e.multiply(recipgamma(a, tmp1, tmp2));
                p = Complex.fromReal(tmp1[0]);
                q = Complex.fromReal(tmp2[0]);
                ka[0] = d.multiply(s).multiply(q).add(a1.multiply(p)).multiply(c);
                f = Complex.newComplex(ka[0]);
                e = Complex.fromReal(a * a);
                p = g.multiply(c).multiply(0.5);
                q = g.invert().multiply(0.5);
                c = Complex.fromReal(1.0);
                d = b.multiply(b);
                ka1[0] = p;
                n = 1;
                do {
                    f = f.multiply(n).add(p).add(q).divide(e.uminus().add(n * n));
                    c = d.multiply(c).divide(n);
                    p.divideEquals(n - a);
                    q.divideEquals(n + a);
                    g = f.multiply(n).uminus().add(p).multiply(c);
                    h = c.multiply(f);
                    ka[0].addEquals(h);
                    ka1[0].addEquals(g);
                    n++;
                } while(h.divide(ka[0]).add(ka1[0].invert().multiply(g.abs())).abs() > 1.0e-15);
                f = ka[0];
                g = ka1[0].divide(b);
            } else {
                Complex expon;
                expon = x.uminus().exp();
                nonexpbesska01(a, x, ka, ka1);
                f = expon.multiply(ka[0]);
                g = expon.multiply(ka1[0]);
            }
            if (rec) {
                x = x.invert().multiply(2.0);
                for (n = 1; n <= na; n++) {
                    h = g.multiply(x).multiply(a + n).add(f);
                    f = g;
                    g = h;
                }
            }
            if (rev) {
                ka1[0] = f;
                ka[0] = g;
            } else {
                ka[0] = f;
                ka1[0] = g;
            }
        }
    }

    public static void nonexpbesska01(double a, double x, double ka[], double ka1[]) {
        if (a == 0.0) {
            nonexpbessk01(x, ka, ka1);
        } else {
            int n, na;
            boolean rec, rev;
            double f, g, h;
            na = 0;
            rev = (a < -0.5);
            if (rev)
                a = -a - 1.0;
            rec = (a >= 0.5);
            if (rec) {
                na = (int) Math.floor(a + 0.5);
                a -= na;
            }
            if (a == -0.5)
                f = g = Math.sqrt(Math.PI / x / 2.0);
            else if (x < 1.0) {
                double expon;
                expon = Math.exp(x);
                besska01(a, x, ka, ka1);
                f = expon * ka[0];
                g = expon * ka1[0];
            } else {
                double b, c, e, p, q;
                c = 0.25 - a * a;
                b = x + x;
                g = 1.0;
                f = 0.0;
                e = Math.cos(a * Math.PI) / Math.PI * x * 1.0e15;
                n = 1;
                do {
                    h = (2.0 * (n + x) * g - (n - 1 + c / n) * f) / (n + 1);
                    f = g;
                    g = h;
                    n++;
                } while (h * n < e);
                p = q = f / g;
                e = b - 2.0;
                do {
                    p = (n - 1 + c / n) / (e + (n + 1) * (2.0 - p));
                    q = p * (1.0 + q);
                    n--;
                } while (n > 0);
                f = Math.sqrt(Math.PI / b) / (1.0 + q);
                g = f * (a + x + 0.5 - p) / x;
            }
            if (rec) {
                x = 2.0 / x;
                for (n = 1; n <= na; n++) {
                    h = f + (a + n) * x * g;
                    f = g;
                    g = h;
                }
            }
            if (rev) {
                ka1[0] = f;
                ka[0] = g;
            } else {
                ka[0] = f;
                ka1[0] = g;
            }
        }
    }

    public static void nonexpbesska01(double a, Complex x, Complex[] ka) {
        Complex[] ka1 = new Complex[1];
        nonexpbesska01(a, x, ka, ka1);
    }

    public static void nonexpbesska01(double a, Complex x, Complex[] ka, Complex[] ka1) {
        if (a == 0.0) {
            nonexpbessk01(x, ka, ka1);
        } else {
            int n, na;
            boolean rec, rev;
            Complex f;
            Complex g;
            Complex h;
            na = 0;
            rev = (a < -0.5);
            if (rev)
                a = -a - 1.0;
            rec = (a >= 0.5);
            if (rec) {
                na = (int) Math.floor(a + 0.5);
                a -= na;
            }
            if (a == -0.5) {
                f = Complex.fromReal(Math.PI).divide(x).multiply(0.5).sqrt();
                g = Complex.newComplex(f);
            }
            else if (x.abs() < 1.0) {
                Complex expon;
                expon = x.exp();
                besska01(a, x, ka, ka1);
                f = expon.multiply(ka[0]);
                g = expon.multiply(ka1[0]);
            } else {
                Complex b;
                double c;
                Complex e;
                Complex p;
                Complex q;
                c = 0.25 - a * a;
                b = x.add(x);
                g = Complex.fromReal(1.0);
                f = Complex.newComplex();
                e = Complex.fromReal(Math.cos(a * Math.PI)).divide(Math.PI).multiply(x).multiply(1.0e15);
                n = 1;
                do {
                    h = x.add(n).multiply(2.0).multiply(g).subtract(f.multiply(n - 1 + c / n)).divide(n + 1);
                    f = g;
                    g = h;
                    n++;
                } while (h.abs() * n < e.abs());
                p = f.divide(g);
                q = Complex.newComplex(p);
                e = b.subtract(2.0);
                do {
                    p = Complex.fromReal(n - 1 + c / n).divide(p.uminus().add(2.0).multiply(n + 1).add(e));
                    q = p.multiply(q.add(1.0));
                    n--;
                } while (n > 0);
                f = b.invert().multiply(Math.PI).sqrt().divide(q.add(1.0));
                g = f.multiply(x.add(a + 0.5).subtract(p)).divide(x);
            }
            if (rec) {
                x = x.invert().multiply(2.0);
                for (n = 1; n <= na; n++) {
                    h = x.multiply(g).multiply(a + n).add(f);
                    f = g;
                    g = h;
                }
            }
            if (rev) {
                ka1[0] = f;
                ka[0] = g;
            } else {
                ka[0] = f;
                ka1[0] = g;
            }
        }
    }

    public static void main(String[] args) {
        double[] ka = new double[1];
        double[] ka1 = new double[1];
        nonexpbesska01(1, 0.5, ka, ka1);
        System.out.println(ka[0]);

        Complex[] kaa = new Complex[1];
        Complex[] kaa1 = new Complex[1];
        besska01(0.1, Complex.newComplex(0.25, 0.25), kaa, kaa1);
        System.out.println(kaa[0]);

        Complex[] kaaa = new Complex[1];
        Complex[] kaaa1 = new Complex[1];
        nonexpbesska01(0.1, Complex.newComplex(0.25, 0.25), kaaa, kaaa1);
        System.out.println(kaaa[0]);

        kaaa = new Complex[1];
        kaaa1 = new Complex[1];
        nonexpbesska01(0.0, Complex.newComplex(5, 5), kaaa, kaaa1);
        System.out.println(kaaa[0]);

        kaaa = new Complex[1];
        kaaa1 = new Complex[1];
        nonexpbesska01(1.0, Complex.newComplex(5, 5), kaaa, kaaa1);
        System.out.println(kaaa[0]);

        kaaa = new Complex[1];
        kaaa1 = new Complex[1];
        nonexpbesska01(0.0, Complex.newComplex(0.5, 0.5), kaaa, kaaa1);
        System.out.println(kaaa[0]);

        kaaa = new Complex[1];
        kaaa1 = new Complex[1];
        nonexpbesska01(0.0, Complex.newComplex(1.0, Math.sqrt(1.25)), kaaa, kaaa1);
        System.out.println(kaaa[0]);

        kaaa = new Complex[1];
        kaaa1 = new Complex[1];
        nonexpbesska01(0.0, Complex.newComplex(1.0, Math.sqrt(1.25) + 0.5), kaaa, kaaa1);
        System.out.println(kaaa[0]);

        kaaa = new Complex[1];
        kaaa1 = new Complex[1];
        nonexpbesska01(0.0, Complex.newComplex(1.0, Math.sqrt(1.25) - 0.5), kaaa, kaaa1);
        System.out.println(kaaa[0]);
    }

}

