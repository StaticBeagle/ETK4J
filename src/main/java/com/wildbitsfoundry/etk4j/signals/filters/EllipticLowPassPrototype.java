package com.wildbitsfoundry.etk4j.signals.filters;

import java.util.ArrayList;
import java.util.List;

import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;

class EllipticLowPassPrototype implements LowPassPrototype {

	@Override
	public
	int getMinOrderNeeded(double fp, double fs, double ap, double as) {
		// Digital Filter Designer's Handbook: With C++ Algorithms by C. Britton
		// Rorabaugh
		double k = fp / fs;
		double kp = Math.sqrt(Math.sqrt(1 - k * k));
		double u = 0.5 * (1 - kp) / (1 + kp);
		double q = u + 2 * Math.pow(u, 5) + 15 * Math.pow(u, 9) + 150 * Math.pow(u, 13);
		double D = (Math.pow(10.0, 0.1 * as) - 1) / (Math.pow(10.0, 0.1 * ap) - 1);

		// Alternative method using elliptic integrals
		// double rt = fp / fs;
		// double kn = Math.sqrt((Math.pow(10.0, 0.1 * ap) - 1) / (Math.pow(10.0, 0.1 *
		// as) - 1));
		// double rtp = Math.sqrt(1 - rt * rt);
		// double knp = Math.sqrt(1 - kn * kn);
		// return compEllipInt1(rt) * compEllipInt1(knp) / (compEllipInt1(rtp) *
		// compEllipInt1(kn));

		double ne = (Math.log10(16.0 * D) / Math.log10(1.0 / q));
		return (int) Math.ceil(ne);
	}

	@Override
	public ZeroPoleGain buildLowPassPrototype(int n, double ap, double as) {
		/*
		 * H. J. Orchard and Alan N. Willson, Jr.,
		 * "Elliptic functions for filter design," IEEE Trans. Circuits Syst., I, vol.
		 * 44, no. 4, pp. 273-287, April 1997.
		 */
		if (n == 1) {
			// filter becomes Chebyshev I
			Complex[] z = new Complex[0];
			Complex[] p = new Complex[1];
			p[0] = new Complex(-Math.sqrt(1.0 / (Math.pow(10.0, ap * 0.1) - 1.0)), 0.0);
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
		Complex j = new Complex(0.0, 1.0);
		Complex mj = j.conj();
		for (int i = 0, m = zeros.length - 1; i < n3; ++i, m = m - 2) {
			double u1 = (2.0 * i + 1.0) * Math.PI / (2.0 * n);
			Complex c = mj.divide(new Complex(-u1, u2).cos());
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
			Complex zero = new Complex(0.0, d / ek.get(0));
			zeros[m] = zero;
			zeros[m - 1] = zero.conj();
		}
		if (n0 == 1) {
			a = 1.0 / Math.sinh(u2);
			for (int en = m1; en >= 1; --en) {
				double k = ek.get(en);
				a = (a - k / a) / (1 + k);
			}
			poles[n - 1] = new Complex(-1.0 / a, 0.0);
		}
		double k = LowPassPrototype.computeGain(zeros, poles);
		if (n % 2 == 0) {
			double eps0 = e.get(0);
			k /= Math.sqrt(1 + eps0 * eps0);
		}
		return new ZeroPoleGain(zeros, poles, k);
	}
}
