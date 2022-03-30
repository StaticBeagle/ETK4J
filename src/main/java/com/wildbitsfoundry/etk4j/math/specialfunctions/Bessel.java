package com.wildbitsfoundry.etk4j.math.specialfunctions;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

import static com.wildbitsfoundry.etk4j.math.specialfunctions.Gamma.recipgamma;

/**
 * Docs https://github.com/JeffBezanson/numal/tree/master/CHAPTER6
 */

public final class Bessel {

    private Bessel() {
    }

    /**
     * Modified Bessel function of the third kind of order zero and one.
     *
     * @param x  Argument at which to evaluate the Bessel function.
     * @param k0 Output: has the value of the modified Bessel function of the third kind of order zero.
     * @param k1 Output: has the value of the modified Bessel function of the third kind of order one.
     */
    public static void bessk01(double x, double[] k0, double[] k1) {
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

    /**
     * Modified Bessel function of the third kind of order zero and one.
     *
     * @param x  Argument at which to evaluate the Bessel function.
     * @param k0 Output: has the value of the Bessel function of order zero.
     * @param k1 Output: has the value of the Bessel function of order one.
     */
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
            d = new Complex(sum0);
            sum1 = d.multiply(-2.0).subtract(1.0);
            c = new Complex(sum1);
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

    /**
     * Exponentially scaled modified Bessel function of the third kind of order zero and one.
     *
     * @param x  Argument at which to evaluate the Bessel function.
     * @param k0 Output: has the value of the Bessel function of order zero.
     * @param k1 Output: has the value of the Bessel function of order one.
     */
    public static void nonexpbessk01(double x, double[] k0, double[] k1) {
        if (x <= 1.5) {
            double expx;
            expx = Math.exp(x);
            bessk01(x, k0, k1);
            k0[0] *= expx;
            k1[0] *= expx;
        } else if (x <= 5.0) {
            int i, r;
            double t2, s1, s2, term1, term2, sqrtexpr, exph2, x2;
            double[] fac = {
                    0.90483741803596, 0.67032004603564,
                    0.40656965974060, 0.20189651799466, 0.82084998623899e-1,
                    0.27323722447293e-1, 0.74465830709243e-2,
                    0.16615572731739e-2, 0.30353913807887e-3,
                    0.45399929762485e-4, 0.55595132416500e-5,
                    0.55739036926944e-6, 0.45753387694459e-7,
                    0.30748798795865e-8, 0.16918979226151e-9,
                    0.76218651945127e-11, 0.28111852987891e-12,
                    0.84890440338729e-14, 0.2098791048793e-15,
                    0.42483542552916e-17
            };
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
            double br, br1, br2, cr, cr1, cr2, ermin1, erplus1, er, f0, f1,
                    expx, y, y2;
            double[] dr = {
                    0.27545e-15, -0.172697e-14,
                    0.1136042e-13, -0.7883236e-13, 0.58081063e-12,
                    -0.457993633e-11, 0.3904375576e-10, -0.36454717921e-9,
                    0.379299645568e-8, -0.450473376411e-7,
                    0.63257510850049e-6, -0.11106685196665e-4,
                    0.26953261276272e-3, -0.11310504646928e-1
            };
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

    /**
     * Exponentially scaled modified Bessel function of the third kind of order zero and one.
     *
     * @param x  Argument at which to evaluate the Bessel function.
     * @param k0 Output: has the value of the Bessel function of order zero.
     * @param k1 Output: has the value of the Bessel function of order one.
     */
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
            double fac[] = {0.90483741803596, 0.67032004603564, 0.40656965974060, 0.20189651799466,
                    0.82084998623899e-1, 0.27323722447293e-1, 0.74465830709243e-2, 0.16615572731739e-2,
                    0.30353913807887e-3, 0.45399929762485e-4, 0.55595132416500e-5, 0.55739036926944e-6,
                    0.45753387694459e-7, 0.30748798795865e-8, 0.16918979226151e-9, 0.76218651945127e-11,
                    0.28111852987891e-12, 0.84890440338729e-14, 0.2098791048793e-15, 0.42483542552916e-17};
            s1 = Complex.fromReal(0.5);
            s2 = new Complex();
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
            double dr[] = {0.27545e-15, -0.172697e-14, 0.1136042e-13, -0.7883236e-13, 0.58081063e-12, -0.457993633e-11,
                    0.3904375576e-10, -0.36454717921e-9, 0.379299645568e-8, -0.450473376411e-7, 0.63257510850049e-6,
                    -0.11106685196665e-4, 0.26953261276272e-3, -0.11310504646928e-1};
            y = x.invert().multiply(10.0).subtract(1.0);
            y2 = y.add(y);
            r = 30;
            br1 = new Complex();
            br2 = new Complex();
            cr1 = new Complex();
            cr2 = new Complex();
            erplus1 = new Complex();
            er = new Complex();
            for (i = 0; i <= 13; i++) {
                r -= 2;
                br = y2.multiply(br1).subtract(br2).add(dr[i]);
                cr = y2.multiply(cr1).subtract(cr2).add(er);
                ermin1 = erplus1.add(r * dr[i]);
                erplus1 = new Complex(er);
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

    /**
     * Modified Bessel function of the third kind of order a and a + 1.
     *
     * @param a   Order of the function.
     * @param x   Argument at which to evaluate the Bessel function.
     * @param ka  Output: has the value of the Bessel function of order a.
     * @param ka1 Output: has the value of the Bessel function of order a + 1.
     */
    public static void besska01(double a, double x, double[] ka, double[] ka1) {
        if (a == 0.0) {
            bessk01(x, ka, ka1);
        } else {
            int n, na = 0;
            boolean rec, rev;
            double f, g, h;
            rev = (a < -0.5);
            if (rev) a = -a - 1.0;
            rec = (a >= 0.5);
            if (rec) {
                na = (int) Math.floor(a + 0.5);
                a -= na;
            }
            if (a == 0.5)
                f = g = Math.sqrt(Math.PI / x / 2.0) * Math.exp(-x);
            else if (x < 1.0) {
                double a1, b, c, d, e, p, q, s;
                double[] tmp1 = new double[1];
                double[] tmp2 = new double[1];
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

    /**
     * Modified Bessel function of the third kind of order a and a + 1.
     *
     * @param a   order of the function.
     * @param x   Argument at which to evaluate the Bessel function.
     * @param ka  Output: has the value of the Bessel function of order a.
     * @param ka1 Output: has the value of the Bessel function of order a + 1.
     */
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
                g = new Complex(f);
            } else if (x.abs() < 1.0) {
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
                f = new Complex(ka[0]);
                e = Complex.fromReal(a * a);
                p = g.multiply(c).multiply(0.5);
                q = g.invert().multiply(0.5);
                c = Complex.fromReal(1.0);
                d = b.multiply(b);
                ka1[0] = new Complex(p);
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
                } while (h.divide(ka[0]).add(ka1[0].invert().multiply(g.abs())).abs() > 1.0e-15);
                f = new Complex(ka[0]);
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

    /**
     * Exponentially scaled modified Bessel function of the third kind of order a and a + 1.
     *
     * @param a   order of the function.
     * @param x   Argument at which to evaluate the Bessel function.
     * @param ka  Output: has the value of the Bessel function of order a.
     * @param ka1 Output: has the value of the Bessel function of order a + 1.
     */
    public static void nonexpbesska01(double a, double x, double[] ka, double[] ka1) {
        if (a == 0.0) {
            nonexpbessk01(x, ka, ka1);
        } else {
            int n, na = 0;
            boolean rec, rev;
            double f, g, h;
            rev = (a < -0.5);
            if (rev) a = -a - 1.0;
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
                f = expon * (ka[0]);
                g = expon * (ka1[0]);
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

    /**
     * Exponentially scaled modified Bessel function of the third kind of order a and a + 1.
     *
     * @param a   order of the function.
     * @param x   Argument at which to evaluate the Bessel function.
     * @param ka  Output: has the value of the Bessel function of order a.
     * @param ka1 Output: has the value of the Bessel function of order a + 1.
     */
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
                g = new Complex(f);
            } else if (x.abs() < 1.0) {
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
                f = new Complex();
                e = Complex.fromReal(Math.cos(a * Math.PI)).divide(Math.PI).multiply(x).multiply(1.0e15);
                n = 1;
                do {
                    h = x.add(n).multiply(2.0).multiply(g).subtract(f.multiply(n - 1 + c / n)).divide(n + 1);
                    f = g;
                    g = h;
                    n++;
                } while (h.abs() * n < e.abs());
                p = f.divide(g);
                q = new Complex(p);
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
}

