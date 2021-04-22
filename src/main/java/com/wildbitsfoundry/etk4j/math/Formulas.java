package com.wildbitsfoundry.etk4j.math;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public class Formulas {
	private Formulas() {}
	
	public static Complex[] quadraticFormula(double a, double b, double c) {
		Complex r1;
		Complex r2;

		if(b == 0 && c == 0) {
			r1 = Complex.newComplex(0.0, 0.0);
			r2 = Complex.newComplex(r1);
			return new Complex[] { r1, r2 };
		}
		
		if(b == 0 && c != 0) {
			double dis = -c / a;
			if(dis < 0) {
				r1 = Complex.newComplex(0.0, Math.sqrt(-dis));
				r2 = r1.conj();
			} else {
				r1 = Complex.newComplex(Math.sqrt(dis), 0.0);
				r2 = Complex.newComplex(-Math.sqrt(dis), 0.0);
			}
			return new Complex[] { r1, r2 };
		}
		
		if(b != 0 && c == 0) {
			r1 = Complex.newComplex(0.0, 0.0);
			r2 = Complex.newComplex(-b / a, 0.0);
			return new Complex[] { r1, r2 };
		}
		
		double dis = b * b - 4 * a * c;

		if (dis < 0) {
			double k = 1 / (2 * a);
			r1 = Complex.newComplex(-b * k, Math.sqrt(-dis) * k);
			r2 = r1.conj();
		} else {
			double q = -0.5 * (b + Math.signum(b) * Math.sqrt(dis));
			r1 = Complex.newComplex(c / q, 0.0);
			r2 = Complex.newComplex(q / a, 0.0);
		}
		return new Complex[] { r1, r2 };
	}
}
