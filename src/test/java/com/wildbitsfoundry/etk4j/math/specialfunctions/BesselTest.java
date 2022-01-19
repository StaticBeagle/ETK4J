package com.wildbitsfoundry.etk4j.math.specialfunctions;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import static com.wildbitsfoundry.etk4j.math.specialfunctions.Bessel.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class BesselTest {

    @Test
    public void testModifiedOfThirdKind() {
        double[] ka = new double[1];
        double[] ka1 = new double[1];
        besska01(1, 0.5, ka, ka1);
        assertEquals(1.6564411200033013, ka[0], 1e-12);
        assertEquals(7.55018355124086, ka1[0], 1e-12);

        besska01(0, 0.5, ka, ka1);
        assertEquals(0.9244190712276656, ka[0], 1e-12);
        assertEquals(1.6564411200033011, ka1[0], 1e-12);

        besska01(0.1, 0.5, ka, ka1);
        assertEquals(0.9300865291314745, ka[0], 1e-12);
        assertEquals(1.8605926626555973, ka1[0], 1e-12);

        besska01(1.1, 0.5, ka, ka1);
        assertEquals(1.8605926626555973, ka[0], 1e-12);
        assertEquals(9.116694244816104, ka1[0], 1e-12);

        besska01(1.0, 0.9, ka, ka1);
        assertEquals(0.7165335787760228, ka[0], 1e-12);
        assertEquals(2.079027149887384, ka1[0], 1e-12);

        besska01(1.1, 1, ka, ka1);
        assertEquals(0.6475743722371923, ka[0], 1e-12);
        assertEquals(1.8472295638769927, ka1[0], 1e-12);

        Complex[] kaComplex = new Complex[1];
        Complex[] ka1Complex = new Complex[1];

        besska01(0.1, new Complex(0.25, 0.25), kaComplex, ka1Complex);
        assertEquals(1.18490903385095, kaComplex[0].real(), 1e-12);
        assertEquals(-0.7275845912845041, kaComplex[0].imag(), 1e-12);
        assertEquals(1.7898874434260907, ka1Complex[0].real(), 1e-12);
        assertEquals(-2.5398443963682755, ka1Complex[0].imag(), 1e-12);

        besska01(1.1, new Complex(0.25, 0.25), kaComplex, ka1Complex);
        assertEquals(1.7898874434260907, kaComplex[0].real(), 1e-12);
        assertEquals(-2.5398443963682755, kaComplex[0].imag(), 1e-12);
        assertEquals(-2.1149015590946636, ka1Complex[0].real(), 1e-12);
        assertEquals(-19.778404686379716, ka1Complex[0].imag(), 1e-12);

        besska01(1.0, new Complex(0.9, 0.5), kaComplex, ka1Complex);
        assertEquals(0.4227757335613111, kaComplex[0].real(), 1e-12);
        assertEquals(-0.4838006788854517, kaComplex[0].imag(), 1e-12);
        assertEquals(0.6092537971963008, ka1Complex[0].real(), 1e-12);
        assertEquals(-1.5246844138859015, ka1Complex[0].imag(), 1e-12);

        besska01(0.1, new Complex(0.5, 0.0), kaComplex, ka1Complex);
        assertEquals(0.9300865291314747, kaComplex[0].real(), 1e-12);
        assertEquals(0.0, kaComplex[0].imag(), 1e-12);
        assertEquals(1.8605926626555973, ka1Complex[0].real(), 1e-12);
        assertEquals(0.0, ka1Complex[0].imag(), 1e-12);

        besska01(0.1, new Complex(5, 5), kaComplex, ka1Complex);
        assertEquals(0.0019472484810153816, kaComplex[0].real(), 1e-12);
        assertEquals(0.002460939365683677, kaComplex[0].imag(), 1e-12);
        assertEquals(0.0022060510043949936, ka1Complex[0].real(), 1e-12);
        assertEquals(0.00249554428463785, ka1Complex[0].imag(), 1e-12);

        besska01(1.1, new Complex(5, 5), kaComplex, ka1Complex);
        assertEquals(0.0022060510043949936, kaComplex[0].real(), 1e-12);
        assertEquals(0.00249554428463785, kaComplex[0].imag(), 1e-12);
        assertEquals(0.002981599444602607, ka1Complex[0].real(), 1e-12);
        assertEquals(0.0025246278873371055, ka1Complex[0].imag(), 1e-12);

        besska01(0.0, new Complex(5, 5), kaComplex, ka1Complex);
        assertEquals(0.0019451630724234032, kaComplex[0].real(), 1e-12);
        assertEquals(0.002460604699964898, kaComplex[0].imag(), 1e-12);
        assertEquals(0.0021595269236999795, ka1Complex[0].real(), 1e-12);
        assertEquals(0.002490313491075621, ka1Complex[0].imag(), 1e-12);

        besska01(1.0, new Complex(5, 5), kaComplex, ka1Complex);
        assertEquals(0.0021595269236999795, kaComplex[0].real(), 1e-12);
        assertEquals(0.002490313491075621, kaComplex[0].imag(), 1e-12);
        assertEquals(0.0028751311553785234, ka1Complex[0].real(), 1e-12);
        assertEquals(0.002526762013440026, ka1Complex[0].imag(), 1e-12);

        besska01(0.0, new Complex(0.5, 0.5), kaComplex, ka1Complex);
        assertEquals(0.5529723109255748, kaComplex[0].real(), 1e-12);
        assertEquals(-0.5996419478565946, kaComplex[0].imag(), 1e-12);
        assertEquals(0.5784533638220992, ka1Complex[0].real(), 1e-12);
        assertEquals(-1.0828582158182143, ka1Complex[0].imag(), 1e-12);

        besska01(1.0, new Complex(0.5, 0.5), kaComplex, ka1Complex);
        assertEquals(0.5784533638220992, kaComplex[0].real(), 1e-12);
        assertEquals(-1.0828582158182143, kaComplex[0].imag(), 1e-12);
        assertEquals(-0.4558373930666556, ka1Complex[0].real(), 1e-12);
        assertEquals(-3.922265107137218, ka1Complex[0].imag(), 1e-12);

        besska01(0.0, new Complex(1.0, Math.sqrt(1.25)), kaComplex, ka1Complex);
        assertEquals(0.027072279535369895, kaComplex[0].real(), 1e-12);
        assertEquals(-0.3561609318675329, kaComplex[0].imag(), 1e-12);
        assertEquals(-0.04212692657375383, ka1Complex[0].real(), 1e-12);
        assertEquals(-0.4390291101780954, ka1Complex[0].imag(), 1e-12);

        besska01(0.0, new Complex(1.0, Math.sqrt(1.25) + 0.5), kaComplex, ka1Complex);
        assertEquals(-0.15802546052735064, kaComplex[0].real(), 1e-12);
        assertEquals(-0.28080861449928934, kaComplex[0].imag(), 1e-12);
        assertEquals(-0.23627880502809914, ka1Complex[0].real(), 1e-12);
        assertEquals(-0.289712004915866, ka1Complex[0].imag(), 1e-12);

        besska01(0.0, new Complex(1.0, Math.sqrt(1.25) - 0.5), kaComplex, ka1Complex);
        assertEquals(0.2576944014432611, kaComplex[0].real(), 1e-12);
        assertEquals(-0.29940972285544365, kaComplex[0].imag(), 1e-12);
        assertEquals(0.28839524876085937, ka1Complex[0].real(), 1e-12);
        assertEquals(-0.44391820627021067, ka1Complex[0].imag(), 1e-12);
    }

    @Test
    public void testModifiedExponentiallyScaledOfThirdKind() {
        double[] ka = new double[1];
        double[] ka1 = new double[1];
        nonexpbesska01(1, 0.5, ka, ka1);
        assertEquals(2.7310097082117863, ka[0], 1e-12);
        assertEquals(12.448148218621055, ka1[0], 1e-12);

        nonexpbesska01(0, 0.5, ka, ka1);
        assertEquals(1.5241093857739092, ka[0], 1e-12);
        assertEquals(2.7310097082117863, ka1[0], 1e-12);

        nonexpbesska01(0.1, 0.5, ka, ka1);
        assertEquals(1.5334534441707164, ka[0], 1e-12);
        assertEquals(3.0675986990288715, ka1[0], 1e-12);

        nonexpbesska01(1.1, 0.5, ka, ka1);
        assertEquals(3.0675986990288715, ka[0], 1e-12);
        assertEquals(15.030887719897754, ka1[0], 1e-12);

        nonexpbesska01(1.0, 0.9, ka, ka1);
        assertEquals(1.7623882196059197, ka[0], 1e-12);
        assertEquals(5.113581646042784, ka1[0], 1e-12);

        nonexpbesska01(1.1, 1, ka, ka1);
        assertEquals(1.7602896486281334, ka[0], 1e-12);
        assertEquals(5.021290556479156, ka1[0], 1e-12);

        Complex[] kaComplex = new Complex[1];
        Complex[] ka1Complex = new Complex[1];

        nonexpbesska01(0.1, new Complex(0.25, 0.25), kaComplex, ka1Complex);
        assertEquals(1.7052889762566512, kaComplex[0].real(), 1e-12);
        assertEquals(-0.5287803645825102, kaComplex[0].imag(), 1e-12);
        assertEquals(3.033653520010336, ka1Complex[0].real(), 1e-12);
        assertEquals(-2.5912423158732216, ka1Complex[0].imag(), 1e-12);

        nonexpbesska01(1.1, new Complex(0.25, 0.25), kaComplex, ka1Complex);
        assertEquals(3.033653520010336, kaComplex[0].real(), 1e-12);
        assertEquals(-2.5912423158732216, kaComplex[0].imag(), 1e-12);
        assertEquals(3.6518982744599535, ka1Complex[0].real(), 1e-12);
        assertEquals(-25.278322042470165, ka1Complex[0].imag(), 1e-12);

        nonexpbesska01(1.0, new Complex(0.9, 0.5), kaComplex, ka1Complex);
        assertEquals(1.4830595396623063, kaComplex[0].real(), 1e-12);
        assertEquals(-0.5457504025035852, kaComplex[0].imag(), 1e-12);
        assertEquals(3.1129798404562656, ka1Complex[0].real(), 1e-12);
        assertEquals(-2.5726086516444813, ka1Complex[0].imag(), 1e-12);

        nonexpbesska01(0.1, new Complex(0.5, 0.0), kaComplex, ka1Complex);
        assertEquals(1.5334534441707168, kaComplex[0].real(), 1e-12);
        assertEquals(0.0, kaComplex[0].imag(), 1e-12);
        assertEquals(3.0675986990288715, ka1Complex[0].real(), 1e-12);
        assertEquals(0.0, ka1Complex[0].imag(), 1e-12);

        nonexpbesska01(0.1, new Complex(5, 5), kaComplex, ka1Complex);
        assertEquals(0.4322110661286333, kaComplex[0].real(), 1e-12);
        assertEquals(-0.17352294380513053, kaComplex[0].imag(), 1e-12);
        assertEquals(0.44803131318721556, ka1Complex[0].real(), 1e-12);
        assertEquals(-0.20889809815211746, ka1Complex[0].imag(), 1e-12);

        nonexpbesska01(1.1, new Complex(5, 5), kaComplex, ka1Complex);
        assertEquals(0.44803131318721556, kaComplex[0].real(), 1e-12);
        assertEquals(-0.20889809815211746, kaComplex[0].imag(), 1e-12);
        assertEquals(0.48482037343635487, ka1Complex[0].real(), 1e-12);
        assertEquals(-0.31804741429978384, ka1Complex[0].imag(), 1e-12);

        nonexpbesska01(0.0, new Complex(5, 5), kaComplex, ka1Complex);
        assertEquals(0.4320756434783003, kaComplex[0].real(), 1e-12);
        assertEquals(-0.17324024390967138, kaComplex[0].imag(), 1e-12);
        assertEquals(0.4453282558120745, ka1Complex[0].real(), 1e-12);
        assertEquals(-0.20249714367635196, ka1Complex[0].imag(), 1e-12);

        nonexpbesska01(1.0, new Complex(5, 5), kaComplex, ka1Complex);
        assertEquals(0.4453282558120745, kaComplex[0].real(), 1e-12);
        assertEquals(-0.2024971436936182, kaComplex[0].imag(), 1e-12);
        assertEquals(0.4806418659019877, ka1Complex[0].real(), 1e-12);
        assertEquals(-0.3028053238162896, ka1Complex[0].imag(), 1e-12);

        nonexpbesska01(0.0, new Complex(0.5, 0.5), kaComplex, ka1Complex);
        assertEquals(1.2740700057330194, kaComplex[0].real(), 1e-12);
        assertEquals(-0.43052443373915755, kaComplex[0].imag(), 1e-12);
        assertEquals(1.6928912856511094, ka1Complex[0].real(), 1e-12);
        assertEquals(-1.1095435340610968, ka1Complex[0].imag(), 1e-12);

        nonexpbesska01(1.0, new Complex(0.5, 0.5), kaComplex, ka1Complex);
        assertEquals(1.6928912856511094, kaComplex[0].real(), 1e-12);
        assertEquals(-1.1095435340610968, kaComplex[0].imag(), 1e-12);
        assertEquals(2.4407655089130444, ka1Complex[0].real(), 1e-12);
        assertEquals(-6.0353940731635705, ka1Complex[0].imag(), 1e-12);

        nonexpbesska01(0.0, new Complex(1.0, Math.sqrt(1.25)), kaComplex, ka1Complex);
        assertEquals(0.9027895696248687, kaComplex[0].real(), 1e-12);
        assertEquals(-0.35734124114045224, kaComplex[0].imag(), 1e-12);
        assertEquals(1.0230661523935811, ka1Complex[0].real(), 1e-12);
        assertEquals(-0.6250311866955051, ka1Complex[0].imag(), 1e-12);

        nonexpbesska01(0.0, new Complex(1.0, Math.sqrt(1.25) + 0.5), kaComplex, ka1Complex);
        assertEquals(0.7827492386732262, kaComplex[0].real(), 1e-12);
        assertEquals(-0.3930346699919637, kaComplex[0].imag(), 1e-12);
        assertEquals(0.8169685720097325, ka1Complex[0].real(), 1e-12);
        assertEquals(-0.6043692159144972, ka1Complex[0].imag(), 1e-12);

        nonexpbesska01(0.0, new Complex(1.0, Math.sqrt(1.25) - 0.5), kaComplex, ka1Complex);
        assertEquals(1.0424992536107143, kaComplex[0].real(), 1e-12);
        assertEquals(-0.2574424664855407, kaComplex[0].imag(), 1e-12);
        assertEquals(1.338125723984699, ka1Complex[0].real(), 1e-12);
        assertEquals(-0.5292382369845088, ka1Complex[0].imag(), 1e-12);
    }
}
