//====================================================================
//
// Project: PureDenoise
// 
// Organization: 
// Florian Luisier
// Biomedical Imaging Group (BIG)
// Ecole Polytechnique Fédérale de Lausanne (EPFL)
// Lausanne, Switzerland
//
// Information:
// http://bigwww.epfl.ch/algorithms/denoise/
//
// Reference:
// [1]  F. Luisier, C. Vonesch, T. Blu, M. Unser, "Fast Interscale Wavelet
//      Denoising of Poisson-corrupted Images", Signal Processing, vol. 90,
//      no. 2, pp. 415-427, February 2010.
// [2]  F. Luisier, "The SURE-LET Approach to Image Denoising," Swiss Federal
//      Institute of Technology Lausanne, EPFL Thesis no. 4566 (2010), 232 p.,
//      January 8, 2010.
// [3]  F. Luisier, C. Vonesch, T. Blu, M. Unser, "Fast Haar-Wavelet Denoising
//      of Multidimensional Fluorescence Microscopy Data", Proceedings of the
//      Sixth IEEE International Symposium on Biomedical Imaging: From Nano to
//      Macro (ISBI'09)}, Boston MA, USA, June 28-July 1, 2009, pp. 310-313.
//
// Conditions of use:
// You'll be free to use this software for research purposes, but you
// should not redistribute it without our consent. In addition, we 
// expect you to include a citation or acknowledgement whenever 
// you present or publish results that are based on it.
//
//====================================================================
package denoisegui;

import additionaluserinterface.WalkBar;
import denoise.Denoising;

public class Timer implements Runnable {
	
	private Thread thread;
	private WalkBar walk;
	private Denoising denoising;
	
	public Timer(WalkBar walk, Denoising denoising) {
		this.walk = walk;
		this.denoising = denoising;
	}
	
	public void start() {
		if (thread == null) {
			thread = new Thread(this);
			thread.setPriority(Thread.MIN_PRIORITY);
			thread.start();
		}
	}
	
	public void stop() {
		walk.progress("End of denoising", 100);
		if (thread != null) {
			thread.interrupt();
			thread = null;
		}
	}
	
	public void run() {
		while(thread == Thread.currentThread()) {
			try {
				Thread.sleep(200);
			}
			catch(InterruptedException e) {
			}
			walk.progress("Denoising in progress", denoising.getProgress());
		}
		
	}
	
}
