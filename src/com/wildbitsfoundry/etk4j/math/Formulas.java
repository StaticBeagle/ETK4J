package com.wildbitsfoundry.etk4j.math;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public class Formulas {
	private Formulas() {}
	
	public static Complex[] quadraticFormula(double a, double b, double c) {
		
		double dis = b * b - 4 * a * c;
		Complex r1;
		Complex r2;
		if (dis < 0) {
			double k = 1 / (2 * a);
			r1 = new Complex(-b * k, Math.sqrt(-dis) * k);
			r2 = r1.conj();
		} else {
			double q = -0.5 * (b + Math.signum(b) * Math.sqrt(dis));
			r1 = new Complex(q / a, 0.0);
			r2 = new Complex(c / q, 0.0);
		}
		return new Complex[] { r1, r2 };
	}
}
