//====================================================================
//
// Project: PureDenoise
// Package: denoise
// Class  : Denoising
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

import java.util.*;
import ij.*;
import ij.plugin.Memory;
import imageware.*;

/*
 * Main class for multidimensional denoising of mixed Poisson-Gaussian noise.
 *
 * @author 	Swiss Federal Institute of Technology Lausanne
 *			Biomedical Imaging Group
 *
 * @version 1.0
 *
 */
public class Denoising {

    private ImageWare input, output;
    private double[] Alpha, Delta, Sigma;
    private int nx, ny, nz;
    private int CYCLESPIN, MULTIFRAME;
    private boolean FRAMEWISE = false;
    private boolean STOP = false;
    private boolean LOG = false;
    private boolean SUCCESS = false;
    private int MAX_THREAD = 20;
    private double MAX_MEMORY, IMAGE_MEMORY;
    private int[] CURRENT = new int[1];
    private Memory memo = new Memory();
    private double progress = 0;	// Status of the progress from 0 to 100.

    /**
     * Constructor of the class Denoising.
     *
     * @param input         input data to be denoised (3D ImageWare object)
     * @param Alpha         double array containing the estimated detector gain 
     *                      for each frame/slice.
     * @param Delta         double array containing the estimated detector of-
     *                      fset for each frame/slice.
     * @param Sigma         double array containing the estimated AWGN standard
     *                      deviation for each frame/slice.
     * @param FRAMEWISE     true->Framewise noise parameters estimation
     *                      false->Global noise parameters estimation.
     * @param CYCLESPIN     number of cycle-spins (CS>0). A high value of CS
     *                      yields a high-quality denoising result, but the
     *                      computation time linearly increases with CS.
     * @param MULTIFRAME    number of adjacent frames/slices to be considered 
     *                      for multi-frame/slices denoising. MF>0 must be odd.
     *
     */
    public Denoising(ImageWare input,
            double[] Alpha,
            double[] Delta,
            double[] Sigma,
            boolean FRAMEWISE,
            int CYCLESPIN,
            int MULTIFRAME) {
        this.input = input;
        this.Alpha = Alpha;
        this.Delta = Delta;
        this.Sigma = Sigma;
        this.FRAMEWISE = FRAMEWISE;
        this.CYCLESPIN = CYCLESPIN;
        this.MULTIFRAME = MULTIFRAME;
        nx = input.getSizeX();
        ny = input.getSizeY();
        nz = input.getSizeZ();
        MAX_MEMORY = memo.maxMemory();//in bytes
        IMAGE_MEMORY = nx * ny * nz * 4.0; //32 bit images in bytes
        MAX_THREAD = (int) Math.max(Math.min(Math.round(0.25*MAX_MEMORY/IMAGE_MEMORY),nz),1);//Limit the number of concurrent threads
        //IJ.log("Max Memory = "+MAX_MEMORY);
        //IJ.log("Image Memory = "+IMAGE_MEMORY);
        //IJ.log("Thread = "+MAX_THREAD);
    }

    /**
     * Method for automatic (either global or framewise) noise parameters
     * estimation.
     */
    final public boolean estimateNoiseParameters() {
        ImageWare buffer = Builder.create(nx, ny, 1, ImageWare.DOUBLE);
        if (FRAMEWISE) {
            CURRENT[0] = 0;
            Vector<ParametersThread> V = new Vector<ParametersThread>();
            V.ensureCapacity(nz);
            ParametersThread pt = null;
            int k = 0, current_prev, k_prev;
            int Z = 0;
            while (Z < nz) {
                Z = Math.min(nz, Z + MAX_THREAD);
                k_prev = k;
                while (CURRENT[0] < Z) {
                    if (STOP) {
                        for (int j = k_prev; j < k; j++) {
                            pt = V.get(j);
                            pt.interrupt();
                        }
                        return false;
                    }
                    current_prev = CURRENT[0];
                    pt = new ParametersThread(input, Alpha, Delta, Sigma, CURRENT, SUCCESS, LOG);
                    pt.start();
                    if (LOG) {
                        IJ.log("Noise parameters estimation for frame no " + (current_prev + 1) + " lauched.");
                    }
                    pt.setPriority(Thread.MIN_PRIORITY);
                    V.insertElementAt(pt, k);
                    while (current_prev == CURRENT[0]) {
                        IJ.wait(1);
                    }
                    k++;
                }
                for (int i = k_prev; i < k; i++) {
                    pt = V.get(i);
                    if (STOP) {
                        for (int j = k_prev; j < k; j++) {
                            pt = V.get(j);
                            pt.interrupt();
                        }
                        return false;
                    }
                    try {
                        pt.join();
                        if (!pt.SUCCESS) {
                            for (int j = k_prev; j < k; j++) {
                                pt = V.get(j);
                                pt.interrupt();
                            }
                            return false;
                        }
                    } catch (Exception ex) {
                        IJ.log("Thread no " + i + ": " + pt.getState());
                        IJ.log("Error: thread cannot be ended.");
                    }
                    double p = 100 * (i + 1.0) / nz;
                    IJ.showStatus("Individual noise parameters estimation: " + IJ.d2s(p, 2) + "%");
                    progress = p;
                }
            }
            IJ.showStatus("");
            if (LOG) {
                IJ.log("---------------------------------------------------");
            }
        } else {
            IJ.showStatus("Global noise parameters estimation...");
            progress = 0;
            int c0 = (nz - 1) / 2;
            int zs = (int) Math.max(c0 - 1, 0);
            int ze = (int) Math.min(c0 + 1, nz - 1);
            buffer = Operations.averageSubStack(input, zs, ze);
            progress = 10;
            double[] noiseParams = Operations.estimateNoiseParams(buffer, 4);
            progress = 90;
            Alpha[0] = (ze - zs + 1) * noiseParams[0];
            Delta[0] = noiseParams[1];
            Sigma[0] = Math.sqrt(ze - zs + 1) * noiseParams[2];
            if (Alpha[0]<=0 || Sigma[0] < 0) {
                return false;
            }
            for (int i = 0; i < nz; i++) {
                Alpha[i] = Alpha[0];
                Delta[i] = Delta[0];
                Sigma[i] = Sigma[0];
            }
            progress = 100;
            IJ.showStatus("");
        }

        return true;
    }

    /**
     * Method for parallel denoising of multidimensional images.
     */
    final public void perform() {
        progress = 0;
        int[] offset = new int[2];
        int[] iters = Operations.numOfDyadicIters(nx, ny);
        double[][][] array;
        double progress0 = 2;
        ImageWare buffer = Builder.create(nx, ny, 1, ImageWare.DOUBLE);
        output = Builder.create(nx, ny, nz, ImageWare.DOUBLE);
        for (int i = 0; i < nz; i++) {
            input.getXY(0, 0, i, buffer);
            buffer.subtract(Delta[i]);
            buffer.divide(Alpha[i]);
            input.putXY(0, 0, i, buffer);
        }
        progress = progress+progress0;
        java.util.Random rand = new java.util.Random(0);
        for (int cs = 0; cs < CYCLESPIN; cs++) {
            IJ.showStatus("Denoising in progress for cycle-spin no " + (cs + 1) + ": " + IJ.d2s(0, 2) + "%");
            if (LOG) {
                IJ.log("-------------- Cycle-spin no " + (cs + 1) + " --------------");
            }
            offset[0] = (int) (cs > 0 ? Math.floor(rand.nextDouble() * nx) : 0);
            offset[1] = (int) (cs > 0 ? Math.floor(rand.nextDouble() * ny) : 0);
            //------------------------------------------------------------------
            ImageWare out = Builder.create(nx, ny, nz, ImageWare.DOUBLE);
            CURRENT[0] = 0;
            Vector<DenoisingThread> V = new Vector<DenoisingThread>();
            V.ensureCapacity(nz);
            DenoisingThread dt = null;
            int k = 0, current_prev, k_prev;
            int Z = 0;
            while (Z < nz) {
                Z = Math.min(nz, Z + MAX_THREAD);
                k_prev = k;
                while (CURRENT[0] < Z) {
                    if (STOP) {
                        output = null;
                        for (int j = k_prev; j < k; j++) {
                            dt = V.get(j);
                            dt.interrupt();
                        }
                        return;
                    }
                    current_prev = CURRENT[0];
                    int[] index = Operations.getAdjacentIndex(current_prev, nz, MULTIFRAME);
                    double[] bufferSigma = Operations.restrictArray(Operations.divide(Sigma, Alpha), index[0], index[1]);
                    ImageWare Coef = Builder.create(nx, ny, index[1] - index[0] + 1, ImageWare.DOUBLE);
                    input.getXYZ(0, 0, index[0], Coef);
                    array = new double[nx][ny][bufferSigma.length];
                    Coef.getBlockXYZ(offset[0], offset[1], 0, array, ImageWare.PERIODIC);
                    Coef = Builder.create(array, ImageWare.DOUBLE);
                    dt = new DenoisingThread(Coef, out, Alpha[current_prev], Delta[current_prev], bufferSigma, iters, index[2], CURRENT, LOG);
                    dt.start();
                    if (LOG) {
                        IJ.log("Denoising of frame no " + (current_prev + 1) + " launched.");
                    }
                    dt.setPriority(Thread.MIN_PRIORITY);
                    V.insertElementAt(dt, k);
                    while (current_prev == CURRENT[0]) {
                        IJ.wait(1);
                    }
                    k++;
                }
                for (int i = k_prev; i < k; i++) {
                    dt = (DenoisingThread) V.get(i);
                    if (STOP) {
                        output = null;
                        for (int j = k_prev; j < k; j++) {
                            dt = V.get(j);
                            dt.interrupt();
                        }
                        return;
                    }
                    try {
                        dt.join();
                    } catch (Exception ex) {
                        IJ.log("Thread no " + i + ": " + dt.getState());
                        IJ.log("Error: thread cannot be ended.");
                    }
                    double p = 100*(i+1.0)/nz;
                    IJ.showStatus("Denoising in progress for cycle-spin no " + (cs + 1) + ": " + IJ.d2s(p, 2) + "%");
                    progress = progress0+(cs+(i + 1.0)/nz)/CYCLESPIN*(100-2*progress0);
                }
            }
            array = new double[nx][ny][nz];
            out.getBlockXYZ(-offset[0], -offset[1], 0, array, ImageWare.PERIODIC);
            output.add(Builder.create(array, ImageWare.DOUBLE));
        }
        output.divide(CYCLESPIN);
        progress = progress+progress0;
    }

    /**
     * Method to set the type of noise parmeters estimation (global or framewise).
     *
     * @param framewise  type of noise parmeters estimation:
     *                   true->Framewise.
     *                   false->Global.
     */
    final public void setFramewise(boolean framewise) {
        FRAMEWISE = framewise;
    }

    /**
     * Method to set the number of cycle-spin(s) for reducing potential blocking
     * artifacts.
     *
     * @param cyclespin  number of cycle-spins (cyclespin>0).
     */
    final public void setCycleSpins(int cyclespin) {
        CYCLESPIN = cyclespin;
    }

    /**
     * Method to interrupt the current denoising task.
     *
     * @param stop  true->Stop the current denoising task.
     */
    final public void setStop(boolean stop) {
        STOP = stop;
    }

    /**
     * Method to display some messages related to the current denoising task.
     *
     * @param log    true->Enable display.
     *               false->Disable display.
     */
    final public void setLog(boolean log) {
        LOG = log;
    }

    /**
     * Method to set the number of adjacent frames/slices for multi-frame/slices
     * denoising.
     *
     * @param multiframe  number of adjacent frames/slices for multi-frame/slices
     *                    denoising.
     */
    final public void setMultiFrame(int multiframe) {
        MULTIFRAME = multiframe;
    }

    /**
     * Method to set the value of the detector gain for each frame/slice.
     *
     * @param alpha  array containing the value of the detector gain for each frame/slice.
     */
    final public void setAlpha(double[] alpha) {
        for (int i = 0; i < nz; i++) {
            Alpha[i] = alpha[i];
        }
    }

    /**
     * Method to set the value of the detector offset for each frame/slice.
     *
     * @param delta  array containing the value of the detector offset for each frame/slice.
     */
    final public void setDelta(double[] delta) {
        for (int i = 0; i < nz; i++) {
            Delta[i] = delta[i];
        }
    }

    /**
     * Method to set the standard deviation of the additive-white-Gaussian-noise
     * (AWGN) for each frame/slice.
     *
     * @param sigma  array containing the standard deviation of the AWGN for each
     *               frame/slice.
     */
    final public void setSigma(double[] sigma) {
        for (int i = 0; i < nz; i++) {
            Sigma[i] = sigma[i];
        }
    }

    /**
     * Method to get the value of the detector gain for each frame/slice.
     *
     * @return  Alpha  array containing the value of the detector gain for each frame/slice.
     */
    final public double[] getAlpha() {
        return Alpha;
    }

    /**
     * Method to get the value of the detector offset for each frame/slice.
     *
     * @return Delta  array containing the value of the detector offset for each frame/slice.
     */
    final public double[] getDelta() {
        return Delta;
    }

    /**
     * Method to get the standard deviation of the additive-white-Gaussian-noise
     * (AWGN) for each frame/slice.
     *
     * @return Sigma  array containing the standard deviation of the AWGN for each
     *                frame/slice.
     */
    final public double[] getSigma() {
        return Sigma;
    }

    /**
     * Method to get the output of the denoising task.
     *
     * @return output  3D ImageWare containing the denoised image.
     */
    final public ImageWare getOutput() {
        return output;
    }

    /**
     * Method to get the maximum number of parallel threads launched.
     *
     * @return MAX_THREAD  the maximum number of parallel threads launched.
     */
    final public int getMaxThread() {
        return MAX_THREAD;
    }

    /**
     * Method to get the level of progression during the denoising process.
     *
     * @return progress level between 0 and 100.
     */
    final public double getProgress() {
        return progress;
    }

    /**
     * Method to detect whether the denoising process has been stopped by the user.
     *
     * @return STOP boolean.
     */
    final public boolean getStop() {
        return STOP;
    }
}