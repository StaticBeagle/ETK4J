package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.complex.Complex;

interface LowPassPrototype {

	int getMinOrderNeeded(double fp, double fs, double ap, double as);

	ZeroPoleGain buildLowPassPrototype(int n, double ap, double as);

	default double getScalingFrequency(double wp, double ws) {
		return wp;
	}

	static double computeGain(Complex[] zeros, Complex[] poles) {
		// Compute gain k
		Complex knum = new Complex(1.0, 0.0);
		for (Complex zero : zeros) {
			knum.multiplyEquals(zero.uminus());
		}
		Complex kden = new Complex(1.0, 0.0);
		for (Complex pole : poles) {
			kden.multiplyEquals(pole.uminus());
		}
		kden.divideEquals(knum);
		return kden.real();
	}
}
