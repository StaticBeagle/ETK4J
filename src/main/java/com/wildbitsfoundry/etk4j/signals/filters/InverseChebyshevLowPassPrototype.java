package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;

class InverseChebyshevLowPassPrototype implements LowPassPrototype {

	@Override
	public
	int getMinOrderNeeded(double fp, double fs, double ap, double as) {
		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		double amax = Math.pow(10, ap * 0.1) - 1;
		double amin = Math.pow(10, as * 0.1) - 1;

		double ne = MathETK.acosh(Math.sqrt(amin / amax)) / MathETK.acosh(ws / wp);
		return (int) Math.ceil(ne);
	}

	@Override
	public
	ZeroPoleGain buildLowPassPrototype(int n, double ap, double as) {
		double eps = 1.0 / Math.sqrt(Math.pow(10, as * 0.1) - 1);

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
				poles[i].divideEquals(Math.pow(poles[i].abs(), 2));
			}
		} else {
			for (int k = -(n - 1) >> 1, i = 0; k <= (n - 1) >> 1; ++k, ++i) {
				double phik = 180.0 * k * nInv;
				poles[i] = new Complex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
				poles[i].divideEquals(Math.pow(poles[i].abs(), 2));
			}
		}

		Complex[] zeros = new Complex[n % 2 == 0 ? n : n - 1];
		for (int k = 0; k < zeros.length;) {
			Complex zero = new Complex(0.0, -1.0 / Math.cos(0.5 * Math.PI * (k + 1) * nInv));
			zeros[k++] = zero;
			zeros[k++] = zero.conj();
		}
		double k = LowPassPrototype.computeGain(zeros, poles);
		return new ZeroPoleGain(zeros, poles, k);
	}

	@Override
	public
	double getScalingFrequency(double wp, double ws) {
		return ws;
	}

}
