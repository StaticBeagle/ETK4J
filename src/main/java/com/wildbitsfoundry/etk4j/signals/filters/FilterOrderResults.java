package com.wildbitsfoundry.etk4j.signals.filters;

public final class FilterOrderResults {

    private FilterOrderResults() {
    }

    public static class OrderAndCutoffFrequency {
        private int n;
        private double wn;

        OrderAndCutoffFrequency(int n, double wn) {
            this.n = n;
            this.wn = wn;
        }

        public int getOrder() {
            return n;
        }

        public double getCutoffFrequency() {
            return wn;
        }
    }

    public static class OrderAndCutoffFrequencies {
        private int n;
        private double wn0;
        private double wn1;

        OrderAndCutoffFrequencies(int n, double wn0, double wn1) {
            this.n = n;
            this.wn0 = wn0;
            this.wn1 = wn1;
        }

        public int getOrder() {
            return n;
        }

        public double getLowerCutoffFrequency() {
            return wn0;
        }

        public double getUpperCutoffFrequency() {
            return wn1;
        }
    }
}
