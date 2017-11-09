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

	public static Complex[] cubicFormula(double a, double b, double c, double d) {

		double a2 = a * a;
		double b2 = b * b;
		double c2 = c * c;
		double abc = a * b * c;
		double disc = 18.0 * abc * d - 4.0 * b2 * b * d + b2 * c2 - 4.0 * a * c2 * c - 27.0 * a2 * d * d;
		double disc0 = b2 - 3.0 * a * c;
		double disc1 = 2.0 * b2 * b - 9.0 * abc + 27.0 * a2 * d;

		Complex k1 = new Complex(-0.5, 0.866025403784439);
		Complex k2 = new Complex(-0.5, -0.866025403784439);
		Complex r = new Complex(-27.0 * a2 * disc, 0.0).sqrt();
		Complex C = r.add(disc1).divide(2d).pow(1.0 / 3.0);
		double k = -1.0 / (3.0 * a);
		return new Complex[] { C.add(b).add(C.pow(-1).multiply(disc0)).multiply(k),
				C.multiply(k1).add(b).add(C.multiply(k1).pow(-1).multiply(disc0)).multiply(k),
				C.multiply(k2).add(b).add(C.multiply(k2).pow(-1).multiply(disc0)).multiply(k) };
	}
}
