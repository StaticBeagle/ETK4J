package com.wildbitsfoundry.etk4j.signals.filters;

import java.util.ArrayList;
import java.util.List;

import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;

public enum ApproximationType {
    CHEBYSHEV {
        @Override
        int getMinOrderNeeded(double fp, double fs, double ap, double as) {
            double wp = 2 * Math.PI * fp;
            double ws = 2 * Math.PI * fs;
            double amax = Math.pow(10, ap * 0.1) - 1;
            double amin = Math.pow(10, as * 0.1) - 1;
            double ne = MathETK.acosh(Math.sqrt(amin / amax)) / MathETK.acosh(ws / wp);
            return (int) Math.ceil(ne);
        }

        @Override
        ZeroPoleGain buildLowPassPrototype(int n, double ap, double as) {
            double eps = Math.sqrt(Math.pow(10, ap * 0.1) - 1);

            double a = 1.0 / n * MathETK.asinh(1 / eps);
            double sinha = Math.sinh(a);
            double cosha = Math.cosh(a);

            Complex[] poles = new Complex[n];
            final double pid = Math.PI / 180.0;
            final double nInv = 1.0 / n;
            if (n % 2 == 0) {
                for (int k = (-n >> 1) + 1, i = 0; k <= n >> 1; ++k, ++i) {
                    double phik = nInv * (180.0 * k - 90.0);
                    poles[i] = Complex.newComplex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
                }
            } else {
                for (int k = -(n - 1) >> 1, i = 0; k <= (n - 1) >> 1; ++k, ++i) {
                    double phik = 180.0 * k * nInv;
                    poles[i] = Complex.newComplex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
                }
            }
            Complex[] zeros = new Complex[0];
            double k = RationalFunction.calculateGain(zeros, poles);
            if (n % 2 == 0) {
                k /= Math.sqrt(1.0 + eps * eps);
            }
            return new ZeroPoleGain(zeros, poles, k);
        }
    },

    BUTTERWORTH {
        @Override
        int getMinOrderNeeded(double fp, double fs, double ap, double as) {
            double wp = 2 * Math.PI * fp;
            double ws = 2 * Math.PI * fs;
            double amax = Math.pow(10, ap * 0.1) - 1;
            double amin = Math.pow(10, as * 0.1) - 1;

            double ne = Math.log10(amin / amax) / (2 * Math.log10(ws / wp));
            return (int) Math.ceil(ne);
        }

        @Override
        ZeroPoleGain buildLowPassPrototype(int n, double ap, double as) {
            double eps = Math.sqrt(Math.pow(10, ap * 0.1) - 1);
            double wb = Math.pow(eps, -1.0 / n);

            final double pid = Math.PI / 180.0;
            final double nInv = 1.0 / n;
            Complex[] poles = new Complex[n];
            if (n % 2 == 0) {
                for (int k = (-n >> 1) + 1, i = 0; k <= n >> 1; ++k, ++i) {
                    double phik = nInv * (180.0 * k - 90.0);
                    poles[i] = Complex.newComplex(-Math.cos(phik * pid), Math.sin(phik * pid));
                    poles[i].multiplyEquals(wb);
                }
            } else {
                for (int k = -(n - 1) >> 1, i = 0; k <= (n - 1) >> 1; ++k, ++i) {
                    double phik = nInv * 180.0 * k;
                    poles[i] = Complex.newComplex(-Math.cos(phik * pid), Math.sin(phik * pid));
                    poles[i].multiplyEquals(wb);
                }
            }

            Complex[] zeros = new Complex[0];
            double k = RationalFunction.calculateGain(zeros, poles);
            return new ZeroPoleGain(zeros, poles, k);
        }
    },

    INVERSE_CHEBYSHEV {
        @Override
        int getMinOrderNeeded(double fp, double fs, double ap, double as) {
            double wp = 2 * Math.PI * fp;
            double ws = 2 * Math.PI * fs;
            double amax = Math.pow(10, ap * 0.1) - 1;
            double amin = Math.pow(10, as * 0.1) - 1;

            double ne = MathETK.acosh(Math.sqrt(amin / amax)) / MathETK.acosh(ws / wp);
            return (int) Math.ceil(ne);
        }

        @Override
        ZeroPoleGain buildLowPassPrototype(int n, double ap, double as) {
            double eps = 1.0 / Math.sqrt(Math.pow(10, as * 0.1) - 1);

            double aPass = Math.pow(10, (0.1 * Math.abs(ap)));
            double aStop = Math.pow(10, (0.1 * Math.abs(as)));
            double wc = Math.cosh(1.0 / n * MathETK.acosh(Math.sqrt((aStop - 1.0) / (aPass - 1.0))));

            double a = 1.0 / n * MathETK.asinh(1 / eps);
            double sinha = Math.sinh(a);
            double cosha = Math.cosh(a);

            Complex[] poles = new Complex[n];
            final double pid = Math.PI / 180.0;
            final double nInv = 1.0 / n;
            if (n % 2 == 0) {
                for (int k = (-n >> 1) + 1, i = 0; k <= n >> 1; ++k, ++i) {
                    double phik = nInv * (180.0 * k - 90.0);
                    poles[i] = Complex.newComplex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
                    poles[i].divideEquals(Math.pow(poles[i].abs(), 2));
                    poles[i].multiplyEquals(wc);
                }
            } else {
                for (int k = -(n - 1) >> 1, i = 0; k <= (n - 1) >> 1; ++k, ++i) {
                    double phik = 180.0 * k * nInv;
                    poles[i] = Complex.newComplex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
                    poles[i].divideEquals(Math.pow(poles[i].abs(), 2));
                    poles[i].multiplyEquals(wc);
                }
            }

            Complex[] zeros = new Complex[n % 2 == 0 ? n : n - 1];
            for (int k = 0; k < zeros.length; ) {
                Complex zero = Complex.fromImaginary(-1.0 / Math.cos(0.5 * Math.PI * (k + 1) * nInv));
                zero.multiplyEquals(wc);
                zeros[k++] = zero;
                zeros[k++] = zero.conj();

            }
            double k = RationalFunction.calculateGain(zeros, poles);
            return new ZeroPoleGain(zeros, poles, k);
        }
    },
    ELLIPTIC {
        @Override
        int getMinOrderNeeded(double fp, double fs, double ap, double as) {
            // 	Digital Filter Designer's Handbook: With C++ Algorithms by C. Britton Rorabaugh
            double k = fp / fs;
            double kp = Math.sqrt(Math.sqrt(1 - k * k));
            double u = 0.5 * (1 - kp) / (1 + kp);
            double q = u + 2 * Math.pow(u, 5) + 15 * Math.pow(u, 9) + 150 * Math.pow(u, 13);
            double D = (Math.pow(10.0, 0.1 * as) - 1) / (Math.pow(10.0, 0.1 * ap) - 1);

            //  Alternative method using elliptic integrals
            //	double rt = fp / fs;
            //	double kn = Math.sqrt((Math.pow(10.0, 0.1 * ap) - 1) / (Math.pow(10.0, 0.1 * as) - 1));
            //	double rtp = Math.sqrt(1 - rt * rt);
            //	double knp = Math.sqrt(1 - kn * kn);
            //	return compEllipInt1(rt) * compEllipInt1(knp) / (compEllipInt1(rtp) * compEllipInt1(kn));

            double ne = (Math.log10(16.0 * D) / Math.log10(1.0 / q));
            return (int) Math.ceil(ne);
        }

        @Override
        ZeroPoleGain buildLowPassPrototype(int n, double ap, double as) {
            if (n == 1) {
                // filter becomes Chebyshev I
                Complex[] z = new Complex[0];
                Complex[] p = new Complex[1];
                p[0] = Complex.fromReal(-Math.sqrt(1.0 / (Math.pow(10.0, ap * 0.1) - 1.0)));
                double k = -p[0].real();
                return new ZeroPoleGain(z, p, k);
            }

            double dbn = Math.log(10.0) * 0.05;
            int n0 = (int) MathETK.rem(n, 2);
            int n3 = (n - n0) >> 1;
            double apn = dbn * ap;
            double asn = dbn * as;

            List<Double> e = new ArrayList<>();
            e.add(Math.sqrt(2.0 * Math.exp(apn) * Math.sinh(apn)));

            List<Double> g = new ArrayList<>();
            g.add(e.get(0) / Math.sqrt(Math.exp(2 * asn) - 1));

            double v = g.get(0);
            int m2 = 0;
            while (v > 1.0e-150) {
                v = (v / (1.0 + Math.sqrt(1 - v * v)));
                v *= v;
                ++m2;
                g.add(v);
            }

            int m1 = 0;
            List<Double> ek = new ArrayList<>(m1);
            for (int i = 0; i < 10; ++i) {
                m1 = m2 + i;
                while (ek.size() <= m1) {
                    ek.add(0.0);
                }
                ek.set(m1, 4.0 * Math.pow((g.get(m2) / 4.0), Math.pow(2.0, i) / n));
                if (ek.get(m1) < 1.0e-14) {
                    break;
                }
            }

            for (int en = m1; en >= 1; --en) {
                ek.set(en - 1, 2.0 * Math.sqrt(ek.get(en)) / (1.0 + ek.get(en)));
            }

            double a = 0.0;
            for (int en = 1; en <= m2; ++en) {
                a = (1.0 + g.get(en)) * e.get(en - 1) * 0.5;
                e.add(a + Math.sqrt(a * a + g.get(en)));
            }

            double u2 = Math.log((1 + Math.sqrt(1 + Math.pow(e.get(m2), 2))) / e.get(m2)) / n;
            Complex[] zeros = new Complex[n % 2 != 0 ? n - 1 : n];
            Complex[] poles = new Complex[n];
            Complex j = Complex.fromImaginary(1.0);
            Complex mj = j.conj();
            for (int i = 0, m = zeros.length - 1; i < n3; ++i, m = m - 2) {
                double u1 = (2.0 * i + 1.0) * Math.PI / (2.0 * n);
                Complex c = mj.divide(Complex.newComplex(-u1, u2).cos());
                double d = 1.0 / Math.cos(u1);
                for (int en = m1; en >= 1; --en) {
                    double k = ek.get(en);
                    c = c.subtract(c.invert().multiply(k));
                    c.divideEquals(1 + k);
                    d = (d + k / d) / (1 + k);
                }
                Complex pole = c.invert();
                poles[m] = pole;
                poles[m - 1] = pole.conj();
                Complex zero = Complex.fromImaginary(d / ek.get(0));
                zeros[m] = zero;
                zeros[m - 1] = zero.conj();
            }
            if (n0 == 1) {
                a = 1.0 / Math.sinh(u2);
                for (int en = m1; en >= 1; --en) {
                    double k = ek.get(en);
                    a = (a - k / a) / (1 + k);
                }
                poles[n - 1] = Complex.fromReal(-1.0 / a);
            }
            double k = RationalFunction.calculateGain(zeros, poles);
            if (n % 2 == 0) {
                double eps0 = e.get(0);
                k /= Math.sqrt(1 + eps0 * eps0);
            }
            return new ZeroPoleGain(zeros, poles, k);
        }
    };

    abstract int getMinOrderNeeded(double fp, double fs, double ap, double as);

    abstract ZeroPoleGain buildLowPassPrototype(int n, double ap, double as);
}