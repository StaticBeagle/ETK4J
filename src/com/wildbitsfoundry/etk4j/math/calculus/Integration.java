package com.wildbitsfoundry.etk4j.math.calculus;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.util.NumArrays;

public final class Integration {
	private Integration() {}
	
	public static double[] cummulativeTrapz(double[] a) {
		final int length = a.length;
		double[] result = new double[length];
		for(int i = 1; i < length; ++i) {
			result[i] = result[i - 1] + (a[i] + a[i - 1]) * 0.5;
		}
		return result;
	}
	
	public static double trapz(double[] a) {
		final int length = a.length;
		double result = 0.0;
		for(int i = 1; i < length; ++i) {
			result += (a[i] + a[i - 1]) * 0.5;
		}
		return result;
	}
	
	public static void main(String[] args) {
		double[] a = new double[] {1, 2, 3, 4 ,5 ,6};
		
		System.out.println(Arrays.toString(cummulativeTrapz(a)));
		System.out.println((trapz(a)));
		System.out.println(Arrays.toString(NumArrays.cummulativeSum(a)));
		System.out.println(NumArrays.sum(a));
		System.out.println(NumArrays.rms(a));
		
	}
}
