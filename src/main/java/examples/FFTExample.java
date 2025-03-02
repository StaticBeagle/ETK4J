package examples;

import java.util.ArrayList;
import java.util.List;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.signal.fft.FFT;

public class FFTExample {
	
    private static double[] real = { 0.559677760695923, 0.811961988413724, -0.489763299068234, -0.778067179542341,
            -0.706735964715945, 0.789019721173435, -0.700001555401106, 0.700489895029539, -0.718488129867960,
            0.574342509973089, 0.685516617700664, 0.817398984407037, -0.267529687111015, -0.224422049737993,
            -0.606118364173603, -0.757421636418101, -0.554515146424656, -0.440790757929746, -0.409781781881146,
            0.0141205685354671, 0.280704586278388, -0.452294658802835, -0.681799182902633, 0.864458500245397,
            -0.435253596623234, -0.0249062599816210, 0.177507536424897, 0.918235637008785, 0.659079582021293,
            -0.0474811926324961, -0.285224755579650, 0.863443049844534 };

    private static double[] imag = { -0.627553075609633, 0.546655823856343, -0.946010653090666, 0.243006473934334,
            -0.257873950960501, 0.253985941117167, 0.549619389302165, 0.278368441048426, -0.896404939976776,
            -0.588221017267259, -0.988610592161680, -0.667865653695379, 0.555211953225667, 0.657655783982322,
            -0.248329940445516, 0.462132700569790, 0.937594247545697, 0.863758857076813, 0.381492588056045,
            -0.151998214979134, -0.929646800306250, 0.0416082390322683, -0.215594991256026, -0.976012577178314,
            0.847837107109384, -0.0561434192544805, 0.531367878780680, -0.296828765390418, -0.0170543373258241,
            -0.222993692923477, 0.775753049706212, -0.142807967492327 };

    public static void main(String[] args) {

    	exampleArrayOfDoubles();
    	System.out.println();
    	exampleWithArrayOfComplex();
    }
    
    private static void exampleArrayOfDoubles() {
        int n = 32;
        
        System.out.println("Two arrays representinng the real and complex part respectively");
    	System.out.println("//////////////////////////////////////////////////////////////////");
        System.out.println("Input data\t\tFFT()\t\t\tIFFT()");
        System.out.println("//////////////////////////////////////////////////////////////////");
        List<StringBuilder> table = new ArrayList<>();

        for (int i = 0; i < n; ++i) {
        	StringBuilder sb = new StringBuilder();
        	sb.append(String.format("(%.4f, %.4f)", real[i], imag[i]));
        	table.add(sb);
        }

        // Perform FFT
        FFT fft = new FFT(n);
        fft.direct(real, imag);
        for (int i = 0; i < n; ++i) {
        	table.get(i).append(String.format("\t(%.4f, %.4f)", real[i], imag[i]));
        }

        // Perform inverse FFT
        fft.inverse(real, imag);
        for (int i = 0; i < n; ++i) {
        	table.get(i).append(String.format("\t(%.4f, %.4f)", real[i], imag[i]));
        }
        
        table.forEach(System.out::println);
    }
    
    private static void exampleWithArrayOfComplex() {
    	int n = 32;
    	
        System.out.println("One array of Complex numbers");
    	System.out.println("//////////////////////////////////////////////////////////////////");
        System.out.println("Input data\t\tFFT()\t\t\tIFFT()");
        System.out.println("//////////////////////////////////////////////////////////////////");
        List<StringBuilder> table = new ArrayList<>();
    	
    	Complex[] data = new Complex[n];
    	
    	for(int i = 0; i < n; ++i) {
    		data[i] = new Complex(real[i], imag[i]);
        	StringBuilder sb = new StringBuilder();
        	sb.append(String.format("%s", data[i]));
        	table.add(sb);
    	}
    	
        // Perform FFT
        FFT fft = new FFT(n);
        fft.direct(data);
        for (int i = 0; i < n; ++i) {
        	table.get(i).append(String.format("\t%s", data[i]));
        }

        // Perform inverse FFT
        fft.inverse(data);
        for (int i = 0; i < n; ++i) {
        	table.get(i).append(String.format("\t%s", data[i]));
        }
        
        table.forEach(System.out::println);
    }
}
