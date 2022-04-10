package com.wildbitsfoundry.etk4j.math.specialfunctions;

public final class Elliptic {

    /**
     * Complete Elliptic integral of the first kind ({@code K}).
     * @param k The parameter of the elliptic integral. The absolute value of {@code K(k<sup>2</sup>} has to be less
     *          than 1.
     * @return The complete elliptic integral of the first kind {@code K(k)}.
     * @see <a href="https://www.codeproject.com/Articles/566614/Elliptic-integrals">Elliptic Integrals</a>
     */
    public static double completeEllipticIntegralFirstKind(double k) {
        double sum;
        double term;
        double above;
        double below;
        sum = 1;
        term = 1;

        above = 1;
        below = 2;

        for (int i = 1; i <= 100; i++) {
            term *= above / below;
            sum += Math.pow(k, i) * Math.pow(term, 2);
            above += 2;
            below += 2;
        }
        sum *= 0.5 * Math.PI;
        return sum;
    }

    /**
     * Complete Elliptic integral of the second kind ({@code E}).
     * @param k The parameter of the elliptic integral. The absolute value of {@code E(k<sup>2</sup>} has to be less
     *          than 1.
     * @return The complete elliptic integral of the second kind {@code E(k)}.
     * @see <a href="https://www.codeproject.com/Articles/566614/Elliptic-integrals">Elliptic Integrals</a>
     */
    public static double completeEllipticIntegralSecondKind(double k) {
        double sum;
        double term;
        double above;
        double below;
        sum = 1;
        term = 1;

        above = 1;
        below = 2;

        for (int i = 1; i <= 100; i++) {
            term *= above / below;
            sum -= Math.pow(k, i) * Math.pow(term, 2) / above;
            above += 2;
            below += 2;
        }
        sum *= 0.5 * Math.PI;
        return sum;
    }

    /**
     * Incomplete Elliptic integral of the first kind ({@code K<sub>inc</sub>}).
     * @param angle The angle in rad/s. The angle should be between {@code 0} and {@code &#960;/2}
     * @param k The parameter is taken as {@code k<sup>2</sup>}.
     * @return The incomplete elliptic integral of the first kind {@code K<sub>inc</sub>(k)}.
     * @see <a href="https://www.codeproject.com/Articles/566614/Elliptic-integrals">Elliptic Integrals</a>
     */
    public static double incompleteEllipticIntegralFirstKind(double angle, double k) {
        double result;
        result = Math.sin(angle) * RF(Math.pow(Math.cos(angle), 2), 1 - k * Math.pow(Math.sin(angle), 2), 1);
        return result;
    }

    /**
     * Incomplete Elliptic integral of the second kind ({@code E<sub>inc</sub>}).
     * @param angle The angle in rad/s. The angle should be between {@code 0} and {@code &#960;/2}
     * @param k The parameter is taken as {@code k<sup>2</sup>}.
     * @return The incomplete elliptic integral of the first kind {@code E<sub>inc</sub>(k)}.
     * @see <a href="https://www.codeproject.com/Articles/566614/Elliptic-integrals">Elliptic Integrals</a>
     */
    public static double incompleteEllipticIntegralSecondKind(double angle, double k) {
        double y = 1 - k * Math.pow(Math.sin(angle), 2);
        return Math.sin(angle) * RF(Math.pow(Math.cos(angle), 2), y, 1) + (-1.0 / 3.0) * k * Math.pow(Math.sin(angle),
                (3.0)) * RD(Math.pow(Math.cos(angle), 2), y, 1);
    }

    /**
     * Carlson's Elliptic integral of the first kind
     * @see <a href="http://en.wikipedia.org/wiki/Carlson_symmetric_form#Series_Expansion">Series Expansion</a>
     */
    private static double RF(double x, double y, double z) {
        double dx, dy, dz;
        double lambda;
        double n = 1.0 / 3.0;
        double mean;
        double tmp;
        do {
            lambda = Math.sqrt(x * y) + Math.sqrt(y * z) + Math.sqrt(z * x);
            x = 0.25 * (x + lambda);
            y = 0.25 * (y + lambda);
            z = 0.25 * (z + lambda);
            mean = (x + y + z) * n;
            tmp = 1 / mean;
            dx = (mean - x) * tmp;
            dy = (mean - y) * tmp;
            dz = (mean - z) * tmp;
        } while (Math.max(Math.max(Math.abs(dx), Math.abs(dy)), Math.abs(dz)) > 1e-7);
        double e2 = dx * dy - dz * dz;
        double e3 = dx * dy * dz;
        double c1 = 1.0 / 24.0;
        double c2 = 0.1;
        double c3 = 3.0 / 44.0;
        double c4 = 1.0 / 14.0;

        double result = 1.0 + (c1 * e2 - c2 - c3 * e3) * e2 + c4 * e3;
        return result / Math.sqrt(mean);
    }

    /**
     * Carlson's Elliptic integral of the second kind
     * @see <a href="http://en.wikipedia.org/wiki/Carlson_symmetric_form#Series_Expansion">Series Expansion</a>
     */
    private static double RD(double x, double y, double z) {
        double dx, dy, dz;
        double lambda;
        double mu;
        double muInv;
        double sum = 0.0;
        double pow4 = 1.0;

        do {
            lambda = Math.sqrt(x * y) + Math.sqrt(y * z) + Math.sqrt(z * x);
            sum += pow4 / (Math.sqrt(z) * (z + lambda));

            pow4 *= 0.25;

            x = 0.25 * (x + lambda);
            y = 0.25 * (y + lambda);
            z = 0.25 * (z + lambda);
            mu = (x + y + 3.0 * z) * 0.2;
            muInv = 1.0 / mu;

            dx = 1 - x * muInv;
            dy = 1 - y * muInv;
            dz = 1 - z * muInv;
        } while (Math.max(Math.max(Math.abs(dx), Math.abs(dy)), Math.abs(dz)) > 1e-7);
        double C1 = 3.0 / 14.0;
        double C2 = 1.0 / 6.0;
        double C3 = 9.0 / 22.0;
        double C4 = 3.0 / 26.0;
        double EA = dx * dy;
        double EB = dz * dz;
        double EC = EA - EB;
        double ED = EA - 6.0 * EB;
        double EF = ED + EC + EC;
        double S1 = ED * (-C1 + 0.25 * C3 * ED - 1.50 * C4 * dz * EF);
        double S2 = dz * (C2 * EF + dz * (-C3 * EC + dz * C4 * EA));

        return 3.0 * sum + pow4 * (1.0 + S1 + S2) / (mu * Math.sqrt(mu));
    }
}