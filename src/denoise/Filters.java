package denoise;

/**
 * <hr>
 * <p><b>Plugin of ImageJ:</b><br>
 * Fractional Spline Wavelet + Some Symlets<br>
 *
 * <p><b>Authors:</b><br>
 * Gil Gaillard, Florian Luisier, Daniel Sage, Dimitri Van De Ville
 * <a href="mailto:florian.luisier@a3.epfl.ch?subject=Fractional Spline Wavelet Package">florian.luisier@a3.epfl.ch</a><br>
 * Swiss Federal Institute of Technology Lausanne,
 * Biomedical Imaging Group,
 * CH-1015 Lausanne, Switzerland,
 * <a href="http://bigwww.epfl.ch">http://bigwww.epfl.ch</a><br>
 *
 * <p><b>Original Version:</b><br>
 * April 2003<br>
 *
 * <p><b>Updated Version:</b><br>
 * 2007<br>
 *
 * <p><b>Copyright</b><br>
 * Copyright 2003, Swiss Federal Institute of Technology, Lausanne, Switzerland, (EPFL)<br>
 *
 * <hr>
 *
 * <p><b>Purpose of the class:</b><br>
 * Generate the coefficients of the filter.
 */

final public class Filters extends Object
{
	protected static final int HAAR         = 0;
	protected static final int ORTHONORMAL  = 1;
	protected static final int SYMLET4      = 2;
	protected static final int SYMLET8      = 3;
	protected static final int BSPLINE      = 4;
	protected static final int DUAL         = 5;
	protected static final int UNHAAR       = 6;

	private double[] Himag;		// Lowpass  filter - Imaginary part 
	private double[] Hreal;		// Lowpass  filter - Real part 
	private double[] Gimag;		// Highpass filter - Imaginary part 
	private double[] Greal;		// Highpass filter - Real part 
	private double[] A; 		// Autocorrelation A(z)
	private double[] A2; 		// Autocorrelation A(z^2)

	private int     filter = 0;
	private int     size = 0;
	private double  degree = 0.0;
	private double  shift = 0.0;
	private int     GDC = 0;

	/**
	* Prepare the computation the spectrum of the desired filter with the same
	* number of samples as the signal to filter.
	*
	* There are Haar and unnormalized Haar filters, spline filters, and
        * symlet filters (with 4 and 8 vanishing moments), but new filters can be added.
	* After having executed the methods generateSynthesisFilter() and
	* generateAnalysisFilters() the result is in four arrays : 
	* Lowpass  filter in H (Hreal[size], Himag[size])
	* Highpass filter in G (Greal[size], Gimag[size])	 
	*
	* @param 	size   		length of the signal to process
	* @param	filter		HAAR/UNHAAR/ORTHONORMAL/BSPLINE/DUAL/SYMLET4/SYMLET8
	* @param	degree		degree of the splines (fractional degree)
	* @param	shift		shift, for the symmetrical version of the filter shift = 0.0
        * @param        GDC             if GDC=1, compute group-delay compensated filters
	*/
	protected Filters(final int size, final int filter, final double degree, final double shift, final int GDC)
	{
		this.size = size;
		this.filter = filter;
		this.degree = degree;
		this.shift = shift;
		this.GDC = GDC;
		Himag = new double[size];	// H(z) Lowpass filter - Imaginary part
		Hreal = new double[size];	// H(z) Lowpass filter - Real part
		Gimag = new double[size];	// G(z) Highpass filter - Imaginary part
		Greal = new double[size];	// H(z) Highpass filter - Real part
		A     = new double[size];	// A(z)
		A2    = new double[size];	// A(z^2) 
		
		computeAutocorrelation();
		computeAutocorrelation2();
	}

	/**
	* Generate all the analysis filter.
	*/
	protected void generateAnalysisFilters()
	{
		switch(filter) {
		
			case ORTHONORMAL:
				generateBspline();			// H(z)
				normalizeOrthogonal();		// H(z) <- H(z) * sqrt(A(z)/A(z^2))
				generateHighpass();			// G(z) <- 1/z * H(z)
				mirrorHighpass();			// G(z) <- G(-1/z)
				break;
		
			case BSPLINE:
				generateBspline();			// H(z)			
				generateHighpass();			// G(z) <- 1/z * H(z)
				divideAutocorrelation2();	// G(z) <- G(z) / A(z^2)
				mirrorHighpass();			// G(z) <- G(-1/z)
				normalizeBspline();			// H(z) <- H(z) * A(z) / A(z^2)
				break;
		
			case DUAL:
				generateBspline();			// H(z)			
				generateHighpass();			// G(z) <- 1/z * H(z)
				multiplyAutocorrelation();	// G(z) <- G(z) * A(z)
				mirrorHighpass();			// G(z) <- G(-1/z)
				break;
				
			case SYMLET4:
				generateSymlet(4);			// H(z)
				generateHighpass();			// G(z) <- 1/z * H(z)
				mirrorHighpass();			// G(z) <- G(-1/z)
				break;
			
			case SYMLET8:
				generateSymlet(8);			// H(z)
				generateHighpass();			// G(z) <- 1/z * H(z)
				mirrorHighpass();			// G(z) <- G(-1/z)
				break;

			case UNHAAR:
				generateUnhaar();			// H(z), G(z)
				break;				
			
			default:
			    throw new IllegalStateException("Wrong filter name.");
		}
		if(GDC==1){
			generateGDCfilter();
		}		
	}

	/**
	* Generate all the synthesis filter.
	*/
	final protected void generateSynthesisFilters()
	{

		switch(filter) {
		
			case ORTHONORMAL:
				generateBspline();			// H(z)			
				normalizeOrthogonal();		// H(z) <- H(z) * sqrt(A(z)/A(z^2))
				generateHighpass();			// G(z) <- 1/z * H(z)
				mirrorHighpass();			// G(z) <- G(-1/z)
				break;
		
			case BSPLINE:
				generateBspline();			// H(z)			
				generateHighpass();			// G(z) <- 1/z * H(z)
				multiplyAutocorrelation();	// G(z) <- G(z) * A(z)
				mirrorHighpass();			// G(z) <- G(-1/z)
				break;
		
			case DUAL:
				generateBspline();			// H(z)			
				generateHighpass();			// G(z) <- 1/z * H(z)
				divideAutocorrelation2();	// G(z) <- G(z) / A(z^2)
				mirrorHighpass();			// G(z) <- G(-1/z)
				normalizeBspline();			// H(z) <- H(z) * A(z) / A(z^2)
				break;
				
			case SYMLET4:
				generateSymlet(4);			// H(z)
				generateHighpass();			// G(z) <- 1/z * H(z)
				mirrorHighpass();			// G(z) <- G(-1/z)
				break;
			
			case SYMLET8:
				generateSymlet(8);			// H(z)
				generateHighpass();			// G(z) <- 1/z * H(z)
				mirrorHighpass();			// G(z) <- G(-1/z)
				break;
			
			case UNHAAR:
				generateUnhaar();			// H(z)
				for(int n=0;n<size;n++){
					Hreal[n] = Hreal[n]/2.0;
					Himag[n] = Himag[n]/2.0;
					Greal[n] = Greal[n]/2.0;
					Gimag[n] = Gimag[n]/2.0;					
				}
				break;				
			default:
			    throw new IllegalStateException("Wrong filter name.");
		}
	}
	
	/**
	* Generate unnormalized Haar filters.
	*/
	final private void generateUnhaar()
	{
		Hreal[0] = 1.0;
		Himag[0] = 0.0;
		Greal[0] = 1.0;
		Gimag[0] = 0.0;	
		if(size>1){	
			Hreal[size-1] = 1.0;
			Himag[size-1] = 0.0;
			Greal[size-1] = -1.0;
			Gimag[size-1] = 0.0;				
		}
		for(int n=1;n<size-1;n++){
			Hreal[n] = 0.0;	
			Himag[n] = 0.0;
			Greal[n] = 0.0;	
			Gimag[n] = 0.0;					
		}
		FFT1D fft = new FFT1D(size);
		fft.transform(Hreal,Himag,size,0);	
		fft.transform(Greal,Gimag,size,0);	
	}	
	
	/**
	* Generate symlet filters with N vanishing moments.
	*/
	final private void generateSymlet(int N)
	{
		for(int n=0;n<Math.min(2*N,size);n++){
			Hreal[n] = getSymletCoef(n,N);
			Himag[n] = 0.0;	
		}
			
		for(int n=2*N;n<size;n++){
			Hreal[n] = 0.0;	
			Himag[n] = 0.0;	
		}
		FFT1D fft = new FFT1D(size);
		fft.transform(Hreal,Himag,size,0);		
	}
	
	/**
	*/
	final private double getSymletCoef(int n,int N)
	{
		double coef = 0.0;
		switch(N){
			case 4:
				  double[] coef4 = {0.03222310060404,-0.01260396726204,-0.09921954357685,0.29785779560528,0.80373875180592,0.49761866763202,-0.02963552764600,-0.07576571478927};
				  coef = coef4[n];
				  break;
			case 8:
				  double[] coef8 = {0.00188995033276,-0.00030292051472,-0.01495225833705,0.00380875201389,0.04913717967361,-0.02721902991706,-0.05194583810771,0.36444189483533,0.77718575170052,0.48135965125837,-0.06127335906766,-0.14329423835081,0.00760748732492,0.03169508781149,-0.00054213233179,-0.00338241595101};
				  coef = coef8[n];
				  break;				  
		}
		
		return coef;		
	}	

	/**
	* Generate H(z) = B(z, a, T).
	*/
	final private void generateBspline()
	{
		double M = (double)size;		// size of the signal in double
		double b, c;					//temporary variables
		double sqrt2 = Math.sqrt(2.0);
		int halfsize = size/2;
		
		double w, w2;
		double order = degree + 1.0;
		
		for (int k=0; k<size; k++) {
			w2 = Math.PI * (double)k / M;
			w = 2.0 * w2;
			if (k <= halfsize) {
				b = shift * w;
				c = sqrt2 * Math.pow(Math.cos(w2), order);
			}
			else {
				b = shift * (w - 2.0 * Math.PI);
				c = sqrt2 * Math.pow(-Math.cos(w2), order);
			}
			Hreal[k] =  c * Math.cos(b);
			Himag[k] = -c * Math.sin(b);
		                           }
	}

	/**
	* H(z) <- H(z) * sqrt(A(z)/A(z^2)).
	*/
	final private void normalizeOrthogonal()
	{
		double factor;
		for (int k=0;k<size;k++) {
			factor = Math.sqrt(A[k] / A2[k]);
			Hreal[k] *= factor;
			Himag[k] *= factor;
		}
	}

	/**
	* H(z) <- H(z) * A(z) / A(z^2).
	*/
	final private void normalizeBspline()
	{
		double factor;
		for (int k=0;k<size;k++) {
			factor = A[k] / A2[k];
			Hreal[k] *= factor;
			Himag[k] *= factor;
		}
	}

	/**
	* G(z) <- z * H(z).
	*/
	final private void generateHighpass()
	{	
		double M = (double)size;		// size of the signal in double
		double d, e, w;
		for (int k=0;k<size;k++) {
			w = 2.0* Math.PI * (double)k / M;
			d = Math.cos(w);
			e = -Math.sin(w);
			Greal[k] = d * Hreal[k] - e * Himag[k];
			Gimag[k] = d * Himag[k] + e * Hreal[k];
		}
	}

	/**
	* G(z) <- G(-1/z).
	*/
	final private void mirrorHighpass()
	{	
		int halfsize = size/2;
		double swap;
		for (int k=0; k<halfsize; k++) {
			swap = Greal[k+halfsize]; 
			Greal[k+halfsize] =  Greal[k]; 
			Greal[k] =  swap;
			swap = Gimag[k+halfsize]; 
			Gimag[k+halfsize] = -Gimag[k]; 
			Gimag[k] = -swap;
		}
	}

	/**
	* G(z) <- G(z) * A(z).
	*/
	final private void multiplyAutocorrelation()
	{
		for (int k=0; k<size; k++) {
			Greal[k] *= A[k];
			Gimag[k] *= A[k];
		}
	}

	/**
	* G(z) <- G(z) / A(z^2).
	*/
	final private void divideAutocorrelation2()
	{
		for (int k=0; k<size; k++) {
			Greal[k] /= A2[k];
			Gimag[k] /= A2[k];
		}
	}

	/**
	* Calculates the autocorrelation of splines in Fourier domain using
	* Poisson-equivalent expression. 
	*/
	final private void computeAutocorrelation()
	{
		double M = (double)size;
		int N = 100;
		double degree2 =  2.0 * degree + 2.0;
		double nu, sum = 0.0;
		
		for (int k=0; k<size; k++) {
			nu = (double)k / M;
			if (nu == 0) 
				A[k] = 1.0;
			else {
				sum = 0.0;
				for (int l=-N; l<=N; l++) {
					sum += 1.0 / Math.pow( Math.abs(nu + l), degree2);
				                          }
				// Correction terms in O(1/N^(2*alpha+5)
				sum += 2.0 / (degree2-1.0) / Math.pow(N, degree2-1.0);
				sum -= 1.0 / (Math.pow(N, degree2));
				sum += (degree+1.0) * ((1.0/3.0) + 2.0 * nu * nu) / Math.pow(N, degree2+1.0);
				sum -= (degree+1.0) * (degree2+1.0) * nu * nu / Math.pow(N, degree2+2.0);
				// 	Factorization of sin(PI*nu)/PI
				sum *= Math.pow(Math.abs(Math.sin(Math.PI*nu) / Math.PI), degree2);
			 	A[k] = sum;
			      }
		                          }
	}

	/**
	* Compute Autocorrelation A(z^2).
	*/
	private void computeAutocorrelation2()
	{
		int halfsize = size/2;
		for (int k=0; k<halfsize; k++) {
			A2[k] = A[k*2];
			A2[k+halfsize] = A[k*2];
		}
	}
	
	/**
	*/
	private void generateGDCfilter()
	{
		generateW2();		
		Operations.complexMultiplication(Greal,Gimag,Hreal,Himag);	
	}
	
	/**
	*/
	private void generateW2()
	{		
		double w, n0 = 0.0, sqrt2 = Math.sqrt(2.0);;
		int halfsize = size/2;
		switch(filter){
			case HAAR:
				for(int i=0;i<halfsize;i++){
					w = 4*Math.PI*i/size;
					Greal[i] = 0.0;
					Gimag[i] = -sqrt2*Math.sin(w);
					Greal[i+halfsize] = Greal[i];
					Gimag[i+halfsize] = Gimag[i];
				} 				
				break;
			case UNHAAR:
				for(int i=0;i<halfsize;i++){
					w = 4*Math.PI*i/size;
					Greal[i] = 0.0;
					Gimag[i] = -2.0*Math.sin(w);
					Greal[i+halfsize] = Greal[i];
					Gimag[i+halfsize] = Gimag[i];
				} 				
				break;				
			case SYMLET4:
				n0 = -4.0;
				break;
			case SYMLET8:
				n0 = -8.0;
				break;
			default:
				n0 = 0.0;												
		}
		if(filter!=HAAR && filter!=UNHAAR){
			for(int i=0;i<halfsize;i++){
				w = 4*Math.PI*i/size;
				Greal[i] = (Math.cos(n0*w)-Math.cos((n0-1)*w))/sqrt2;
				Gimag[i] = -(Math.sin(n0*w)-Math.sin((n0-1)*w))/sqrt2;
				Greal[i+halfsize] = Greal[i];
				Gimag[i+halfsize] = Gimag[i];				
			} 			
		}
	}		

	/**
	* Return the real parts of the Highpass filter.	 
	*/
	protected double [] getRealHighpassFilter()
	{
		return Greal;
	}

	/**
	* Return the imaginary parts of the Highpass filter.	 
	*/
	protected double [] getImaginaryHighpassFilter()
	{
		return Gimag;
	}

	/**
	* Return the real parts of the Lowpass filter.	 
	*/
	protected double [] getRealLowpassFilter()
	{
		return Hreal;
	}

	/**
	* Return the imaginary parts of the Lowpass filter.
	*/
	protected double [] getImaginaryLowpassFilter()
	{
		return Himag;
	}

} // end of class
