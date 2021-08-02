package com.wildbitsfoundry.etk4j.statistics.regression;

interface RegressionModel {
	double R2();
	double SSE();
	double SSR();
	double SST();
	double[] beta();
	double[] residuals();
	double normOfResiduals();
}
