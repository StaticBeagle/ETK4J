package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;

class ChebyshevLowPassPrototype implements LowPassPrototype {

	@Override
	public int getMinOrderNeeded(double fp, double fs, double ap, double as) {
		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		double amax = Math.pow(10, ap * 0.1) - 1;
		double amin = Math.pow(10, as * 0.1) - 1;
		double ne = MathETK.acosh(Math.sqrt(amin / amax)) / MathETK.acosh(ws / wp);
		return (int) Math.ceil(ne);
	}

	@Override
	public ZeroPoleGain buildLowPassPrototype(int n, double ap, double as) {
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
				poles[i] = new Complex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
			}
		} else {
			for (int k = -(n - 1) >> 1, i = 0; k <= (n - 1) >> 1; ++k, ++i) {
				double phik = 180.0 * k * nInv;
				poles[i] = new Complex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
			}
		}
		Complex[] zeros = new Complex[0];
		double k = LowPassPrototype.computeGain(zeros, poles);
		if (n % 2 == 0) {
			k /= Math.sqrt(1.0 + eps * eps);
		}
		return new ZeroPoleGain(zeros, poles, k);
	}

}
