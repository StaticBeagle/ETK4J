package com.wildbitsfoundry.etk4j.math;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public class Formulas {
	private Formulas() {}
	
	public static Complex[] quadraticFormula(double a, double b, double c) {
		Complex r1;
		Complex r2;

		double dis = b * b - 4 * a * c;
		double inv = 1.0 / (2.0 * a);
		if(dis < 0) {
			r1 = new Complex(-b * inv + 0.0 ,  Math.sqrt(-dis) * inv);
			r2 = r1.conj();
		} else {
			r1 = Complex.fromReal((-b + Math.sqrt(dis)) * inv + 0.0);
			r2 = Complex.fromReal((-b - Math.sqrt(dis)) * inv + 0.0);
		}
		return new Complex[] { r1, r2 };
	}
}
