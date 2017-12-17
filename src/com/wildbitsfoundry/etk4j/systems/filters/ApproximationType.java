package com.wildbitsfoundry.etk4j.systems.filters;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.systems.TransferFunction;
import com.wildbitsfoundry.etk4j.systems.filters.AnalogFilter.LowPassPrototype;
import com.wildbitsfoundry.etk4j.util.NumArrays;

public enum ApproximationType {
	CHEBYSHEV {
		@Override
		int getMinOrderNeeded(double fp, double fs, double ap, double as) {
			double wp = 2 * Math.PI * fp;
			double ws = 2 * Math.PI * fs;
			double amax = Math.pow(10, ap * 0.1) - 1;
			double amin = Math.pow(10, as * 0.1) - 1;

			double L = MathETK.acosh(Math.sqrt(amin / amax)) / MathETK.acosh(ws / wp);

			return (int) Math.ceil(L);
		}

		@Override
		LowPassPrototype buildLowPassPrototype(int n, double ap) {
			double eps = Math.sqrt(Math.pow(10, ap * 0.1) - 1);

			double a = 1.0 / n * MathETK.asinh(1 / eps);
			double sinha = Math.sinh(a);
			double cosha = Math.cosh(a);

			final double pid = Math.PI / 180.0;
			Complex[] poles = new Complex[n];
			if (n % 2 == 0) {
				int i = 0;
				for (double k : NumArrays.linsteps(-n * 0.5 + 1.0, n * 0.5, 1)) {
					double phik = 180.0 * (k / n) - 90.0 / n;
					poles[i++] = new Complex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
				}
			} else {
				int i = 0;
				for (double k : NumArrays.linsteps(-(n - 1) * 0.5, (n - 1) * 0.5, 1)) {
					double phik = 180.0 * (k / n);
					poles[i++] = new Complex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
				}
			}

			double N = 1;
			for (int k = 0; k < (int) Math.ceil(n / 2.0); ++k) {
				if (poles[k].imag() != 0.0) {
					N *= poles[k].real() * poles[k].real() + poles[k].imag() * poles[k].imag();
				} else {
					N *= -poles[k].real();
				}
			}
			TransferFunction tf = new TransferFunction(N, poles);
			return new LowPassPrototype(eps, tf);
		}

		@Override
		double getBandPassAp(double ap, double as1, double as2) {
			return ap;
		}

		@Override
		double getBandPassBW(int n, double eps, double Q, double w0, double omega) {
			return Q;
		}

		@Override
		double getBandStopAp(double amax, double amin) {
			return amax;
		}

		@Override
		double getBandStopBW(int n, double eps, double Q, double w0, double omegas) {
			return Q;
		}
	},
	
	BUTTERWORTH {
		@Override
		int getMinOrderNeeded(double fp, double fs, double ap, double as) {
			double wp = 2 * Math.PI * fp;
			double ws = 2 * Math.PI * fs;
			double amax = Math.pow(10, ap * 0.1) - 1;
			double amin = Math.pow(10, as * 0.1) - 1;

			double L = Math.log10(amin / amax) / (2 * Math.log10(ws / wp));

			return (int) Math.ceil(L);
		}

		@Override
		LowPassPrototype buildLowPassPrototype(int n, double ap) {
			double eps = Math.sqrt(Math.pow(10, ap * 0.1) - 1);

			final double pid = Math.PI / 180.0;
			Complex[] poles = new Complex[n];
			if (n % 2 == 0) {
				int i = 0;
				for (double k : NumArrays.linsteps(-n * 0.5 + 1.0, n * 0.5, 1)) {
					double phik = 180.0 * (k / n) - 90.0 / n;
					poles[i++] = new Complex(-Math.cos(phik * pid), Math.sin(phik * pid));
				}
			} else {
				int i = 0;
				for (double k : NumArrays.linsteps(-(n - 1) * 0.5, (n - 1) * 0.5, 1)) {
					double phik = 180.0 * (k / n);
					poles[i++] = new Complex(-Math.cos(phik * pid), Math.sin(phik * pid));
				}
			}
			TransferFunction tf = new TransferFunction(new Complex[0], poles);
			return new LowPassPrototype(eps, tf);
		}

		@Override
		double getBandPassAp(double ap, double as1, double as2) {
			return ap;
		}

		@Override
		double getBandPassBW(int n, double eps, double Q, double w0, double omega) {
			return Q / Math.pow(eps, -1.0 / n) / w0;
		}

		@Override
		double getBandStopAp(double amax, double amin) {
			return amax;
		}

		@Override
		double getBandStopBW(int n, double eps, double Q, double w0, double omegas) {
			return Q * Math.pow(eps, -1.0 / n) / w0;
		}
	},
	
	INVERSE_CHEBYSHEV
	{
		@Override
		int getMinOrderNeeded(double fp, double fs, double ap, double as) {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		LowPassPrototype buildLowPassPrototype(int n, double ap) {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		double getBandPassAp(double ap, double as1, double as2) {
			return as2;
		}

		@Override
		double getBandPassBW(int n, double eps, double Q, double w0, double omega) {
			return Q / omega;
		}

		@Override
		double getBandStopAp(double amax, double amin) {
			return amin;
		}

		@Override
		double getBandStopBW(int n, double eps, double Q, double w0, double omegas) {
			return Q * omegas;
		}
	};
	abstract int getMinOrderNeeded(double fp, double fs, double ap, double as);
	abstract LowPassPrototype buildLowPassPrototype(int n, double ap);
	abstract double getBandPassAp(double ap, double as1, double as2);
	abstract double getBandPassBW(int n, double eps, double Q, double w0, double omega);
	abstract double getBandStopAp(double amax, double amin);
	abstract double getBandStopBW(int n, double eps, double Q, double w0, double omegas);
}
