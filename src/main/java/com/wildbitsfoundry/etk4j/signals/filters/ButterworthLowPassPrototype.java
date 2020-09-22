package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.complex.Complex;

class ButterworthLowPassPrototype implements LowPassPrototype {

	@Override
	public int getMinOrderNeeded(double fp, double fs, double ap, double as) {
		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		double amax = Math.pow(10, ap * 0.1) - 1;
		double amin = Math.pow(10, as * 0.1) - 1;

		double ne = Math.log10(amin / amax) / (2 * Math.log10(ws / wp));
		return (int) Math.ceil(ne);
	}

	@Override
	public ZeroPoleGain buildLowPassPrototype(int n, double ap, double as) {
		double eps = Math.sqrt(Math.pow(10, ap * 0.1) - 1);
		double wb = Math.pow(eps, -1.0 / n);

		final double pid = Math.PI / 180.0;
		final double nInv = 1.0 / n;
		Complex[] poles = new Complex[n];
		if (n % 2 == 0) {
			for (int k = (-n >> 1) + 1, i = 0; k <= n >> 1; ++k, ++i) {
				double phik = nInv * (180.0 * k - 90.0);
				poles[i] = new Complex(-Math.cos(phik * pid), Math.sin(phik * pid));
				poles[i].multiplyEquals(wb);
			}
		} else {
			for (int k = -(n - 1) >> 1, i = 0; k <= (n - 1) >> 1; ++k, ++i) {
				double phik = nInv * 180.0 * k;
				poles[i] = new Complex(-Math.cos(phik * pid), Math.sin(phik * pid));
				poles[i].multiplyEquals(wb);
			}
		}

		Complex[] zeros = new Complex[0];
		double k = LowPassPrototype.computeGain(zeros, poles);
		return new ZeroPoleGain(zeros, poles, k);
	}

}
