package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

import static com.wildbitsfoundry.etk4j.util.ComplexArrays.deepCopy;

public class ZeroPoleGain {
	private final Complex[] zeros;
	private final Complex[] poles;
	private final double gain;
	
	public ZeroPoleGain(Complex[] zeros, Complex[] poles, double gain) {
		this.zeros = deepCopy(zeros);
		this.poles = deepCopy(poles);
		this.gain = gain;
	}

	public Complex[] getZeros() {
		return deepCopy(zeros);
	}

	public Complex[] getPoles() {
		return deepCopy(poles);
	}

	public double getGain() {
		return gain;
	}
}
