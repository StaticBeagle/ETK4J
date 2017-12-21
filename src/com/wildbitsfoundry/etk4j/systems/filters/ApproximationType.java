package com.wildbitsfoundry.etk4j.systems.filters;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.systems.TransferFunction;
import com.wildbitsfoundry.etk4j.systems.filters.AnalogFilter.LowPassPrototype;
import com.wildbitsfoundry.etk4j.util.NumArrays;
import com.wildbitsfoundry.etk4j.util.Tuples.Tuple2;

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
				for (double k : NumArrays.linsteps(-n * 0.5 + 1.0, 1, n * 0.5)) {
					double phik = 180.0 * (k / n) - 90.0 / n;
					poles[i++] = new Complex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
				}
			} else {
				int i = 0;
				for (double k : NumArrays.linsteps(-(n - 1) * 0.5, 1, (n - 1) * 0.5)) {
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

		@Override
		double getLowPassGainFactor(int n, double eps, double wp) {
			return wp;
		}

		@Override
		double getLowPassAttenuation(double ap, double as) {
			return ap;
		}

		@Override
		double getHighPassAttenuation(double ap, double as) {
			return ap;
		}

		@Override
		double getHighPassGainFactor(int n, double eps, double wp, double ws) {
			return wp;
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
				for (double k : NumArrays.linsteps(-n * 0.5 + 1.0, 1, n * 0.5)) {
					double phik = 180.0 * (k / n) - 90.0 / n;
					poles[i++] = new Complex(-Math.cos(phik * pid), Math.sin(phik * pid));
				}
			} else {
				int i = 0;
				for (double k : NumArrays.linsteps(-(n - 1) * 0.5, 1, (n - 1) * 0.5)) {
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

		@Override
		double getLowPassGainFactor(int n, double eps, double wp) {
			return Math.pow(eps, -1.0 / n) * wp;
		}

		@Override
		double getLowPassAttenuation(double ap, double as) {
			return ap;
		}

		@Override
		double getHighPassAttenuation(double ap, double as) {
			return ap;
		}

		@Override
		double getHighPassGainFactor(int n, double eps, double wp, double ws) {
			return wp / Math.pow(eps, -1.0 / n);
		}
	},
	
	INVERSE_CHEBYSHEV
	{
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
			double eps = 1.0 / Math.sqrt(Math.pow(10, ap * 0.1) - 1);

			double a = 1.0 / n * MathETK.asinh(1 / eps);
			double sinha = Math.sinh(a);
			double cosha = Math.cosh(a);

			Complex[] poles = new Complex[n];
			List<Complex> zeros = new ArrayList<>(n);
			if (n % 2 == 0) {
				for (int k = -n >> 1 + 1, i = 0; k <= n >> 1; ++k, ++i) {
					double phik = 180.0 * (1.0 * k / n) - 90.0 / n;
					Tuple2<Complex, Optional<Complex>> pz = this.calcPZ(phik, cosha, sinha, k, n);
					poles[i] = pz.Item1;
					if(pz.Item2.isPresent()) {
						zeros.add(pz.Item2.get());
					}
				}
			} else {
				for (int k = -(n - 1) >> 1, i = 0; k <= (n - 1) >> 1; ++k, ++i) {
					double phik = 180.0 * (1.0 * k / n);
					Tuple2<Complex, Optional<Complex>> pz = this.calcPZ(phik, cosha, sinha, k, n);
					poles[i] = pz.Item1;
					if(pz.Item2.isPresent()) {
						zeros.add(pz.Item2.get());
					}
				}
			}
			
			double G = 1;
			for (int k = 0; k < (int) Math.ceil(n / 2.0); ++k) {
				if (poles[k].imag() != 0.0) {
					G *= poles[k].real() * poles[k].real() + poles[k].imag() * poles[k].imag();
				} else {
					G *= -poles[k].real();
				}
			}
			double Gp = 1.00;
			for(int k = 0; k <= (int) Math.ceil(zeros.size() >> 1); k = k + 2) {
				Gp *= zeros.get(k).imag();
			}
			Gp *= Gp;
			G /= Gp;
			
			Polynomial num = new Polynomial(zeros.toArray(new Complex[zeros.size()]));
			num.multiplyEquals(G);
			Polynomial den = new Polynomial(poles);
			TransferFunction tf = new TransferFunction(num, den);

			return new LowPassPrototype(eps, tf);
		}

		private Tuple2<Complex, Optional<Complex>> calcPZ(double phik, double cosha, double sinha, int k, int n) {
			final double pid = Math.PI / 180.0;
			Complex pole = new Complex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
			pole.divideEquals(Math.pow(pole.abs(), 2));
			Complex zero = null;
			double phikz = Math.PI / n * (k + 0.5);
			double sign = k < 0 ? -1 : 1;
			if(phikz != 0.5 * Math.PI) {
				zero = new Complex(0.0, sign / Math.cos(phikz));						
			}
			return Tuple2.createTuple(pole, Optional.ofNullable(zero));
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

		@Override
		double getLowPassGainFactor(int n, double eps, double wp) {
			return wp;
		}

		@Override
		double getLowPassAttenuation(double ap, double as) {
			return as;
		}

		@Override
		double getHighPassAttenuation(double ap, double as) {
			return as;
		}

		@Override
		double getHighPassGainFactor(int n, double eps, double wp, double ws) {
			return ws;
		}
	};
	abstract int getMinOrderNeeded(double fp, double fs, double ap, double as);
	abstract LowPassPrototype buildLowPassPrototype(int n, double ap);
	abstract double getBandPassAp(double ap, double as1, double as2);
	abstract double getBandPassBW(int n, double eps, double Q, double w0, double omega);
	abstract double getBandStopAp(double amax, double amin);
	abstract double getBandStopBW(int n, double eps, double Q, double w0, double omegas);
	abstract double getLowPassGainFactor(int n, double eps, double wp);
	abstract double getLowPassAttenuation(double ap, double as);
	abstract double getHighPassAttenuation(double ap, double as);
	abstract double getHighPassGainFactor(int n, double eps, double wp, double ws);
}
