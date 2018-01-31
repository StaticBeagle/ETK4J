package com.wildbitsfoundry.etk4j.math.specialfunctions;

import java.util.ArrayList;
import java.util.List;

public final class Elliptic {
	private Elliptic() {
	}

//	Carlson, B.C. 1977, SIAM Journal on Mathematical Analysis, vol. 8, pp. 231–242. [3]
//			Carlson, B.C. 1987, Mathematics of Computation, vol. 49, pp. 595–606 [4]; 1988, op. cit., vol. 51,
//			pp. 267–280 [5]; 1989, op. cit., vol. 53, pp. 327–333 [6]; 1991, op. cit., vol. 56, pp. 267–280.
//			[7]
//			Bulirsch, R
	// K(k) = RF(0, 1 - k^2), 1);
	
	// Legendre form 
	public static double compEllipInt1(double k) {
		return carlsonsRF(0.0, 1 - k * k, 1.0);
	}

	// K(k) = RF(0, 1 - k^2), 1); Legendre form
	public static double compEllipInt1(double k, double tol) {
		return carlsonsRF(0.0, 1 - k * k, 1.0, tol);
	}

	// Legendre form
	public static double incompEllipInt1(double phi, double k, double tol) {
		double sinPhi = Math.sin(phi);
		double x = Math.pow(Math.cos(phi), 2);
		double y = 1 - k * k * sinPhi * sinPhi;
		return sinPhi * carlsonsRF(x, y, 1.0, tol);
	}

	// Carlson's RF elliptic integral of the first kind
	public static double carlsonsRF(double x, double y, double z) {
		return carlsonsRF(x, y, z, 0.0025);
	}

	// Carlson's RF elliptic integral of the first kind
	public static double carlsonsRF(double x, double y, double z, double tol) {
		double minLim = Double.MIN_VALUE * 5.0; // Five times the min value
		double maxLim = Double.MAX_VALUE * 0.2; // One fifth of max value

		if (Math.min(Math.min(x, y), z) < 0.0 || Math.min(Math.min(x + y, x + z), y + z) < minLim
				|| Math.max(Math.max(x, y), z) > maxLim) {
			throw new IllegalArgumentException("Limits exceeded in carlsonsRF");
		}

		double dx, dy, dz;
		double lambda;
		double n = 1.0 / 3.0;
		double mean;
		double tmp;
		do {
			lambda = Math.sqrt(x * y) + Math.sqrt(y * z) + Math.sqrt(z * x);
			x = 0.25 * (x + lambda);
			y = 0.25 * (y + lambda);
			z = 0.25 * (z + lambda);
			mean = (x + y + z) * n;
			tmp = 1 / mean;
			dx = (mean - x) * tmp;
			dy = (mean - y) * tmp;
			dz = (mean - z) * tmp;
		} while (Math.max(Math.max(Math.abs(dx), Math.abs(dy)), Math.abs(dz)) > tol);
		double e2 = dx * dy - dz * dz;
		double e3 = dx * dy * dz;
		double c1 = 1.0 / 24.0;
		double c2 = 0.1;
		double c3 = 3.0 / 44.0;
		double c4 = 1.0 / 14.0;

		double result = 1.0 + (c1 * e2 - c2 - c3 * e3) * e2 + c4 * e3;
		return result / Math.sqrt(mean);
	}

	// Carlson's elliptic integral of the second kind
	public static double carlsonsRD(double x, double y, double z) {
		return carlsonsRD(x, y, z, 0.0015);
	}

	// Carlson's elliptic integral of the second kind
	public static double carlsonsRD(double x, double y, double z, double tol) {
		double minLim = 2.0 * Math.pow(Double.MAX_VALUE, -2.0 / 3.0);
		double maxLim = 0.1 * tol * Math.pow(Double.MIN_VALUE, -2.0 / 3.0);

		if (Math.min(x, y) < 0.0 || Math.min(x + y, z) < minLim || Math.max(Math.max(x, y), z) > maxLim) {
			throw new IllegalArgumentException("Limits exceeded in carlsonsRD");
		}

		double dx, dy, dz;
		double lambda;
		double mu;
		double muInv;
		double sum = 0.0;
		double pow4 = 1.0;
		
		do {
			lambda = Math.sqrt(x * y) + Math.sqrt(y * z) + Math.sqrt(z * x);
			sum += pow4 / (Math.sqrt(z) * (z + lambda));

			pow4 *= 0.25;

			x = 0.25 * (x + lambda);
			y = 0.25 * (y + lambda);
			z = 0.25 * (z + lambda);
			mu = (x + y + 3.0 * z) * 0.2;
			muInv = 1.0 / mu;

			dx = 1 - x * muInv;
			dy = 1 - y * muInv;
			;
			dz = 1 - z * muInv;
			;
		} while (Math.max(Math.max(Math.abs(dx), Math.abs(dy)), Math.abs(dz)) > tol);
		double C1 = 3.0 / 14.0;
		double C2 = 1.0 / 6.0;
		double C3 = 9.0 / 22.0;
		double C4 = 3.0 / 26.0;
		double EA = dx * dy;
		double EB = dz * dz;
		double EC = EA - EB;
		double ED = EA - 6.0 * EB;
		double EF = ED + EC + EC;
		double S1 = ED * (-C1 + 0.25 * C3 * ED - 1.50 * C4 * dz * EF);
		double S2 = dz * (C2 * EF + dz * (-C3 * EC + dz * C4 * EA));

		return 3.0 * sum + pow4 * (1.0 + S1 + S2) / (mu * Math.sqrt(mu));
	}

	public static double[] landen(double k, double tol) {
		if (k == 0 || k == 1) {
			return new double[] { k };
		}

		List<Double> v = new ArrayList<>();
		while (k > tol) {
			k /= 1 + Math.sqrt(1 - k * k);
			k *= k;
			v.add(k);
		}
		double[] result = new double[v.size()];
		for(int i = 0; i < v.size(); ++i) {
			result[i] = v.get(i).doubleValue();
		}
		return result;
	}

	public static double[] landen(double k, int noIter) {
		if (k == 0 || k == 1) {
			return new double[] { k };
		}

		double[] v = new double[noIter];
		int i = 0;
		while (noIter-- > 0) {
			k /= 1 + Math.sqrt(1 - k * k);
			k *= k;
			v[i++] = k;
		}
		return v;
	}
}
