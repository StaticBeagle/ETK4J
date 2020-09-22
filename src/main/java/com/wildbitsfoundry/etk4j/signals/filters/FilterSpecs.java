package com.wildbitsfoundry.etk4j.signals.filters;

public abstract class FilterSpecs {

	
	public static class LowPassSpecs extends FilterSpecs {

		private double passBandFrequency = 1;
		private double stopBandFrequency = 10;
		private double passBandRipple = 0.2;
		private double stopBandAttenuation = 40;
		
		public double getPassBandFrequency() {
			return passBandFrequency;
		}
		public void setPassBandFrequency(double passBandFrequency) {
			this.passBandFrequency = passBandFrequency;
		}
		public double getStopBandFrequency() {
			return stopBandFrequency;
		}
		public void setStopBandFrequency(double stopBandFrequency) {
			this.stopBandFrequency = stopBandFrequency;
		}
		public double getPassBandRipple() {
			return passBandRipple;
		}
		public void setPassBandRipple(double passBandRipple) {
			this.passBandRipple = passBandRipple;
		}
		public double getStopBandAttenuation() {
			return stopBandAttenuation;
		}
		public void setStopBandAttenuation(double stopBandAttenuation) {
			this.stopBandAttenuation = stopBandAttenuation;
		}
	}
	
	public static class HighPassSpecs extends FilterSpecs {

		private double passBandFrequency = 1;
		private double stopBandFrequency = 10;
		private double passBandRipple = 0.2;
		private double stopBandAttenuation = 40;
		
		public double getPassBandFrequency() {
			return passBandFrequency;
		}
		public void setPassBandFrequency(double passBandFrequency) {
			this.passBandFrequency = passBandFrequency;
		}
		public double getStopBandFrequency() {
			return stopBandFrequency;
		}
		public void setStopBandFrequency(double stopBandFrequency) {
			this.stopBandFrequency = stopBandFrequency;
		}
		public double getPassBandRipple() {
			return passBandRipple;
		}
		public void setPassBandRipple(double passBandRipple) {
			this.passBandRipple = passBandRipple;
		}
		public double getStopBandAttenuation() {
			return stopBandAttenuation;
		}
		public void setStopBandAttenuation(double stopBandAttenuation) {
			this.stopBandAttenuation = stopBandAttenuation;
		}
	}
	
	public static class BandPassSpecs extends FilterSpecs {

		private double lowerPassBandFrequency = 900;
		private double upperPassBandFrequency = 1000;
		private double lowerStopBandFrequency = 90;
		private double upperStopBandFrequency = 10000;
		private double passBandRipple = 0.2;
		private double lowerStopBandAttenuation = 40;
		private double upperStopBandAttenuation = 40;
		
		public double getLowerPassBandFrequency() {
			return lowerPassBandFrequency;
		}
		public void setLowerPassBandFrequency(double lowerPassBandFrequency) {
			this.lowerPassBandFrequency = lowerPassBandFrequency;
		}
		public double getUpperPassBandFrequency() {
			return upperPassBandFrequency;
		}
		public void setUpperPassBandFrequency(double upperPassBandFrequency) {
			this.upperPassBandFrequency = upperPassBandFrequency;
		}
		public double getLowerStopBandFrequency() {
			return lowerStopBandFrequency;
		}
		public void setLowerStopBandFrequency(double lowerStopBandFrequency) {
			this.lowerStopBandFrequency = lowerStopBandFrequency;
		}
		public double getUpperStopBandFrequency() {
			return upperStopBandFrequency;
		}
		public void setUpperStopBandFrequency(double upperStopBandFrequency) {
			this.upperStopBandFrequency = upperStopBandFrequency;
		}
		public double getPassBandRipple() {
			return passBandRipple;
		}
		public void setPassBandRipple(double passBandRipple) {
			this.passBandRipple = passBandRipple;
		}
		public double getLowerStopBandAttenuation() {
			return lowerStopBandAttenuation;
		}
		public void setLowerStopBandAttenuation(double lowerStopBandAttenuation) {
			this.lowerStopBandAttenuation = lowerStopBandAttenuation;
		}
		public double getUpperStopBandAttenuation() {
			return upperStopBandAttenuation;
		}
		public void setUpperStopBandAttenuation(double upperStopBandAttenuation) {
			this.upperStopBandAttenuation = upperStopBandAttenuation;
		}
	}
	
	public static class BandStopSpecs extends FilterSpecs {
		private double lowerPassBandFrequency = 3.6e3;
		private double upperPassBandFrequency = 9.1e3;
		private double lowerStopBandFrequency = 5.45e3;
		private double upperStopBandFrequency = 5.90e3;
		private double passBandRipple = 1.5;
		private double stopBandAttenuation = 40;
		
		public double getLowerPassBandFrequency() {
			return lowerPassBandFrequency;
		}
		public void setLowerPassBandFrequency(double lowerPassBandFrequency) {
			this.lowerPassBandFrequency = lowerPassBandFrequency;
		}
		public double getUpperPassBandFrequency() {
			return upperPassBandFrequency;
		}
		public void setUpperPassBandFrequency(double upperPassBandFrequency) {
			this.upperPassBandFrequency = upperPassBandFrequency;
		}
		public double getLowerStopBandFrequency() {
			return lowerStopBandFrequency;
		}
		public void setLowerStopBandFrequency(double lowerStopBandFrequency) {
			this.lowerStopBandFrequency = lowerStopBandFrequency;
		}
		public double getUpperStopBandFrequency() {
			return upperStopBandFrequency;
		}
		public void setUpperStopBandFrequency(double upperStopBandFrequency) {
			this.upperStopBandFrequency = upperStopBandFrequency;
		}
		public double getPassBandRipple() {
			return passBandRipple;
		}
		public void setPassBandRipple(double passBandRipple) {
			this.passBandRipple = passBandRipple;
		}
		public double getStopBandAttenuation() {
			return stopBandAttenuation;
		}
		public void setStopBandAttenuation(double stopBandAttenuation) {
			this.stopBandAttenuation = stopBandAttenuation;
		}
	}

}
