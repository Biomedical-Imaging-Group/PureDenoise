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
// References:
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
package denoise;

import ij.*;
import imageware.*;

public class ParametersThread extends Thread {

    ImageWare input;
    double[] AlphaHat;
    double[] DeltaHat;
    double[] SigmaHat;
    int[] current;
    int CS = 4;
    boolean SUCCESS;
    boolean LOG;

    protected ParametersThread(ImageWare input,
            double[] AlphaHat,
            double[] DeltaHat,
            double[] SigmaHat,
            int[] current,
            boolean SUCCESS,
            boolean LOG) {
        this.input = input;
        this.AlphaHat = AlphaHat;
        this.DeltaHat = DeltaHat;
        this.SigmaHat = SigmaHat;
        this.current = current;
        this.SUCCESS = SUCCESS;
        this.LOG = LOG;
    }

    @Override
    public void run() {
        int c = current[0];
        current[0] = current[0] + 1;
        ImageWare buffer = Builder.create(input.getWidth(), input.getHeight(), 1, ImageWare.DOUBLE);
        input.getXY(0, 0, c, buffer);
        double[] noiseParams = Operations.estimateNoiseParams(buffer, CS);
        AlphaHat[c] = noiseParams[0];
        DeltaHat[c] = noiseParams[1];
        SigmaHat[c] = noiseParams[2];
        if (LOG) {
            IJ.log("Noise parameters for frame no " + (c + 1) + " estimated.");
        }
        if(AlphaHat[c]<=0 || SigmaHat[c]<0)
            SUCCESS = false;
        else
            SUCCESS = true;
    }
}