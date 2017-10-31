package com.wildbitsfoundry.etk4j.regression;

public interface RegressionModel {
	public double R2();
	public double SSE();
	public double SST();
	public double[] beta();
	public double[] residuals();
}
