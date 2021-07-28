package com.wildbitsfoundry.etk4j.signals.filters;

public abstract class FilterSpecs {



    public static class LowPassSpecs extends FilterSpecs {

        private double passBandFrequency;
        private double stopBandFrequency;
        private double passBandRipple;
        private double stopBandAttenuation;

        public LowPassSpecs() {}

        public LowPassSpecs(LowPassSpecs specs) {
            this.passBandFrequency = specs.passBandFrequency;
            this.stopBandFrequency = specs.stopBandFrequency;
            this.passBandRipple = specs.passBandRipple;
            this.stopBandAttenuation = specs.stopBandAttenuation;
        }

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

        private double passBandFrequency;
        private double stopBandFrequency;
        private double passBandRipple;
        private double stopBandAttenuation;

        public HighPassSpecs() {}

        public HighPassSpecs(HighPassSpecs specs) {
            this.passBandFrequency = specs.passBandFrequency;
            this.stopBandFrequency = specs.stopBandFrequency;
            this.passBandRipple = specs.passBandRipple;
            this.stopBandAttenuation = specs.stopBandAttenuation;
        }

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

        private double lowerPassBandFrequency;
        private double upperPassBandFrequency;
        private double lowerStopBandFrequency;
        private double upperStopBandFrequency;
        private double lowerStopBandAttenuation;
        private double upperStopBandAttenuation;
        private double passBandRipple;
        private double stopBandAttenuation;

        public BandPassSpecs() {}

        public BandPassSpecs(BandPassSpecs specs) {
            this.lowerPassBandFrequency = specs.lowerPassBandFrequency;
            this.upperPassBandFrequency = specs.upperPassBandFrequency;
            this.lowerStopBandFrequency = specs.lowerStopBandFrequency;
            this.upperStopBandFrequency = specs.upperStopBandFrequency;
            this.passBandRipple = specs.passBandRipple;
            this.stopBandAttenuation = specs.stopBandAttenuation;
        }

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

        private double lowerPassBandFrequency;
        private double upperPassBandFrequency;
        private double lowerStopBandFrequency;
        private double upperStopBandFrequency;
        private double passBandRipple;
        private double stopBandAttenuation;

        public BandStopSpecs() {}

        public BandStopSpecs(BandStopSpecs specs) {
            this.lowerPassBandFrequency = specs.lowerPassBandFrequency;
            this.upperPassBandFrequency = specs.upperPassBandFrequency;
            this.lowerStopBandFrequency = specs.lowerStopBandFrequency;
            this.upperStopBandFrequency = specs.upperStopBandFrequency;
            this.passBandRipple = specs.passBandRipple;
            this.stopBandAttenuation = specs.stopBandAttenuation;
        }

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
