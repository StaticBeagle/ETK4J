package com.wildbitsfoundry.etk4j.systems;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public class ZeroPoleGain {
	public Complex[] Zeros;
	public Complex[] Poles;
	public double Gain;
	
	public ZeroPoleGain(Complex[] zeros, Complex[] poles, double gain) {
		this.Zeros = zeros;
		this.Poles = poles;
		this.Gain = gain;
	}
}
