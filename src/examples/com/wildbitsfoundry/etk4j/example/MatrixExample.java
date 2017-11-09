package examples.com.wildbitsfoundry.etk4j.example;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;

public class MatrixExample {
	public static void main(String[] args) {

		double[][] original = { { 1, 2, 3 }, { 0, 4, 5 }, { 1, 0, 6 } };
		Matrix sol = new Matrix(new double[] { 6, 4, 2 }, 3);

		Matrix A = new Matrix(original);
		System.out.printf("A : %n%s%n%n", A);
		
		// Basic matrix operations
		System.out.printf("A.transpose() : %n%s%n%n", A.transpose());
		System.out.printf("A.cofactor() : %n%s%n%n", A.cofactor());
		System.out.printf("A.adjoint() : %n%s%n%n", A.adjoint());
		System.out.printf("A.inv() : %n%s%n%n", A.inv());
		System.out.printf("A.pinv() : %n%s%n%n", A.pinv());
		System.out.printf("A.det() : %n%,4f%n%n", A.det());
		System.out.printf("A.multiply(A) : %n%s%n%n", A.multiply(A));
		System.out.printf("A.add(A) : %n%s%n%n", A.add(A));
		
		// Norms
		System.out.printf("A.norm1() : %n%.4f%n%n", A.norm1());
		System.out.printf("A.norm2() : %n%.4f%n%n", A.norm2());
		System.out.printf("A.normInf() : %n%.4f%n%n", A.normInf());
		System.out.printf("A.normFrob() : %n%.4f%n%n", A.normFrob());
		
		// Sub-matrix operations
		System.out.printf("A.subMatrix(0, 2, 1, 2).transpose() : %n%s%n%n", A.subMatrix(0, 2, 1, 2).transpose());
		System.out.printf("A.subMatrix(1, 2, 0, 2).transpose() : %n%s%n%n", A.subMatrix(1, 2, 0, 2).transpose());
		System.out.printf("A.subMatrix(1, 2, 1, 2) : %n%s%n%n", A.subMatrix(1, 2, 1, 2));
		System.out.printf("A.subMatrix([1, 2], 1, 2) : %n%s%n%n", A.subMatrix(new int[] { 1, 2 }, 1, 2));
		System.out.printf("A.subMatrix(1, 2, [1, 2]) : %n%s%n%n", A.subMatrix(1, 2, new int[] { 1, 2 }));
		System.out.printf("A.subMatrix([1, 2], [1, 2]) : %n%s%n%n",
				A.subMatrix(new int[] { 1, 2 }, new int[] { 1, 2 }));
		System.out.printf("A.subMatrix([0, 2], [1, 2]) : %n%s%n%n",
				A.subMatrix(new int[] { 0, 2 }, new int[] { 1, 2 }));
		
		// Linear system solution
		System.out.printf("A.solve([6, 4, 2]) : %n%s%n%n", A.solve(sol));
		System.out.printf("A(:, 1:2).solve([6, 4, 2]) : %n%s%n%n", A.subMatrix(0, 2, 1, 2).solve(sol));
		System.out.printf("A.eig().getD().diag() : %n%s%n%n", Arrays.toString(A.eig().getD().diag()));

	}
}
