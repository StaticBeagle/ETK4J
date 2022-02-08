package com.wildbitsfoundry.etk4j.signals.filters;

class EllipticOrderCalculationStrategy implements FilterOrderCalculationStrategy{

    @Override
    public double calculateExactOrder(double nat, double gPass, double gStop) {
        // 	Digital Filter Designer's Handbook: With C++ Algorithms by C. Britton Rorabaugh
        double k = 1 / nat;
        double kp = Math.sqrt(Math.sqrt(1 - k * k));
        double u = 0.5 * (1 - kp) / (1 + kp);
        double q = u + 2 * Math.pow(u, 5) + 15 * Math.pow(u, 9) + 150 * Math.pow(u, 13);
        double D = (Math.pow(10.0, 0.1 * gStop) - 1) / (Math.pow(10.0, 0.1 * gPass) - 1);

        //  Alternative method using elliptic integrals
        //	double rt = fp / fs;
        //	double kn = Math.sqrt((Math.pow(10.0, 0.1 * ap) - 1) / (Math.pow(10.0, 0.1 * as) - 1));
        //	double rtp = Math.sqrt(1 - rt * rt);
        //	double knp = Math.sqrt(1 - kn * kn);
        //	return compEllipInt1(rt) * compEllipInt1(knp) / (compEllipInt1(rtp) * compEllipInt1(kn));

        return Math.log10(16.0 * D) / Math.log10(1.0 / q);
    }

    @Override
    public double calculateLowPassWn(int n, LowPassSpecs specs) {
        return specs.getPassBandFrequency();
    }

    @Override
    public double calculateHighPassWn(int n, HighPassSpecs specs) {
        return specs.getPassBandFrequency();
    }

    @Override
    public double[] calculateBandPassWn(int n, BandpassSpecs specs) {
        double[] wn = new double[2];
        wn[0] = specs.getLowerPassBandFrequency();
        wn[1] = specs.getUpperPassBandFrequency();
        return wn;
    }

    @Override
    public double[] calculateBandStopWn(int n, BandStopSpecs specs) {
        double wp1 = specs.getLowerPassBandFrequency();
        double wp2 = specs.getUpperPassBandFrequency();
        return new double[]{wp1, wp2};
    }
}
