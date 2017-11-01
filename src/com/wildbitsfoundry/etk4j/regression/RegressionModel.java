package com.wildbitsfoundry.etk4j.regression;

interface RegressionModel {
	public double R2();
	public double SSE();
	public double SSR();
	public double SST();
	public double[] beta();
	public double[] residuals();
	public double normOfResiduals();
}
