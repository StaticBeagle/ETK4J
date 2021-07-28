package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;
import com.wildbitsfoundry.etk4j.util.Tuples;

import java.util.Arrays;

import static com.wildbitsfoundry.etk4j.math.optimize.minimizers.GoldenSection.goldenSectionMinimizer;

public class Filter {

    protected static FilterOrderResults.OrderAndCutoffFrequency lowPassFilterOrder(LowPassSpecs specs,
                                                                                   FilterOrderCalculationStrategy strategy) {
        double wp = specs.getPassBandFrequency();
        double ws = specs.getStopBandFrequency();
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        double nat = ws / wp;

        double gStop = Math.pow(10, 0.1 * rs);
        double gPass = Math.pow(10, 0.1 * rp);
        int n = strategy.calculateMinOrder(nat, gPass, gStop);
        double wn = strategy.calculateLowPassWn(n, specs);

        return new FilterOrderResults.OrderAndCutoffFrequency(n, wn);
    }

    protected static FilterOrderResults.OrderAndCutoffFrequency highPassFilterOrder(HighPassSpecs specs, FilterOrderCalculationStrategy strategy) {
        double wp = specs.getPassBandFrequency();
        double ws = specs.getStopBandFrequency();
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        double nat = wp / ws;

        double gStop = Math.pow(10, 0.1 * rs);
        double gPass = Math.pow(10, 0.1 * rp);
        int n = strategy.calculateMinOrder(nat, gPass, gStop);
        double wn = strategy.calculateHighPassWn(n, specs);

        return new FilterOrderResults.OrderAndCutoffFrequency(n, wn);
    }

    protected static FilterOrderResults.OrderAndCutoffFrequencies bandPassFilterOrder(FilterSpecs.BandPassSpecs specs,
                                                                                FilterOrderCalculationStrategy strategy) {
        double wp1 = specs.getLowerPassBandFrequency();
        double wp2 = specs.getUpperPassBandFrequency();
        double ws1 = specs.getLowerStopBandFrequency();
        double ws2 = specs.getUpperStopBandFrequency();
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        double w1 = (ws1 * ws1 - wp1 * wp2) / (ws1 * (wp1 - wp2));
        double w2 = (ws2 * ws2 - wp1 * wp2) / (ws2 * (wp1 - wp2));

        double nat = Math.min(Math.abs(w1), Math.abs(w2));

        double gStop = Math.pow(10, 0.1 * rs);
        double gPass = Math.pow(10, 0.1 * rp);

        int n = strategy.calculateMinOrder(nat, gPass, gStop);
        double[] wn = strategy.calculateBandPassWn(n, specs);
        return new FilterOrderResults.OrderAndCutoffFrequencies(n, wn[0], wn[1]);
    }

    protected static FilterOrderResults.OrderAndCutoffFrequencies bandStopFilterOrder(FilterSpecs.BandStopSpecs specs,
                                                                                FilterOrderCalculationStrategy strategy) {
        double wp1 = specs.getLowerPassBandFrequency();
        double wp2 = specs.getUpperPassBandFrequency();
        double ws1 = specs.getLowerStopBandFrequency();
        double ws2 = specs.getUpperStopBandFrequency();
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        double[] wp = new double[2];
        // maximize the pass band
        // https://github.com/scipy/scipy/blob/master/scipy/signal/filter_design.py
        wp[0] = goldenSectionMinimizer(Filter::bandStopObjMinimize, wp1, ws1 - 1e-12,
                1e-5, 500, specs, strategy, 0);

        wp[1] = goldenSectionMinimizer(Filter::bandStopObjMinimize, ws2 + 1e-12, wp2,
                1e-5, 500, specs, strategy, 1);

        wp1 = wp[0];
        wp2 = wp[1];
        double w1 = ws1 * (wp[0] - wp[1]) / (ws1 * ws1 - wp[0] * wp[1]);
        double w2 = ws2 * (wp[0] - wp[1]) / (ws2 * ws2 - wp[0] * wp[1]);

        double nat = Math.min(Math.abs(w1), Math.abs(w2));

        double gStop = Math.pow(10, 0.1 * rs);
        double gPass = Math.pow(10, 0.1 * rp);
        int n = strategy.calculateMinOrder(nat, gPass, gStop);

        specs.setLowerPassBandFrequency(wp1);
        specs.setUpperPassBandFrequency(wp2);
        double[] wn = strategy.calculateBandStopWn(n, specs);

        Arrays.sort(wn);
        return new FilterOrderResults.OrderAndCutoffFrequencies(n, wn[0], wn[1]);
    }

    private static double bandStopObjMinimize(double wp, Object... params) {
        BandStopSpecs specs = (BandStopSpecs) params[0];
        FilterOrderCalculationStrategy type = (FilterOrderCalculationStrategy) params[1];
        Integer index = (Integer) params[2];

        double[] passb = new double[]{specs.getLowerPassBandFrequency(), specs.getUpperPassBandFrequency()};
        double[] stopb = new double[]{specs.getLowerStopBandFrequency(), specs.getUpperStopBandFrequency()};
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        passb[index] = wp;

        double w1 = (stopb[0] * (passb[0] - passb[1]) /
                (stopb[0] * stopb[0] - passb[0] * passb[1]));
        double w2 = (stopb[1] * (passb[0] - passb[1]) /
                (stopb[1] * stopb[1] - passb[0] * passb[1]));

        double gStop = Math.pow(10, 0.1 * rs);
        double gPass = Math.pow(10, 0.1 * rp);
        double nat = Math.min(Math.abs(w1), Math.abs(w2));
        return type.calculateExactOrder(nat, gPass, gStop);
    }
}
