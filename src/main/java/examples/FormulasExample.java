package examples;

import static com.wildbitsfoundry.etk4j.math.Formulas.quadraticFormula;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public class FormulasExample {
	
	public static void main(String[] args) {
		double a = 1.0;
		double b = 3.0;
		double c = -4;
		
		// a * x^2 + b * x + c 
		// roots should be [-4, 1]
		System.out.println("Real distinct roots");
		for(Complex root : quadraticFormula(a, b, c)) {
			System.out.println(root);
		}
		
		// imaginary roots 
		// roots should be [-0.5000 + 0.8660i, -0.5000 - 0.8660i]
		System.out.printf("%nComplex roots %n");
		for(Complex root : quadraticFormula(1, 1, 1)) {
			System.out.println(root);
		}
		
	}
}
