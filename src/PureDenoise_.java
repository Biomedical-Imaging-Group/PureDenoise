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
// expect you to include a citation or acknowledgment whenever 
// you present or publish results that are based on it.
//
//====================================================================

import ij.IJ;
import ij.ImagePlus;
import ij.LookUpTable;
import ij.Macro;
import ij.WindowManager;
import ij.process.ImageProcessor;
import imageware.Builder;
import imageware.Display;
import imageware.ImageWare;

import java.awt.image.ColorModel;
import java.awt.image.IndexColorModel;
import java.util.StringTokenizer;

import denoise.Denoising;
import denoise.Operations;
import denoisegui.DenoiseDialog;

public class PureDenoise_ {

    private double Alpha = 1.0, Delta = 0.0, Sigma = 0.0;
    private double Imax = 0.0, Imin = 0.0;
    private String titleIn = "";
    private String titleOut = "";
    private double[] AlphaHat, DeltaHat, SigmaHat;	// estimated gain & offset of the detectors, AGWN std
    private Denoising denoising;
    private ImageWare original = null;		// Input image
    private ImageWare output = null;		// Output image
    private int CS = 1;
    private int NBFRAME = 0;
    private int Nmin = 16;
    private boolean FRAMEWISE = false;
    private boolean COLOR = false;
    private boolean LOG = false;
    private byte[] red, green, blue;
    private int nx = 0;
    private int ny = 0;
    private int nz = 0;
    private int nxe = 0;
    private int nye = 0;
    private int[] Ext = new int[2];

    public PureDenoise_() {
        if (IJ.versionLessThan("1.25")) {
            return;
        }
        String options = Macro.getOptions();
        if (options != null) {
            double time = System.currentTimeMillis();
            ImagePlus impSource = WindowManager.getCurrentImage();
            if (impSource == null) {
                return;
            }

            int type = impSource.getType();
            titleIn = impSource.getTitle();
            titleOut = "Denoised-" + titleIn;
            LookUpTable lut = impSource.createLut();
            red = lut.getReds();
            green = lut.getGreens();
            blue = lut.getBlues();

            nx = impSource.getWidth();
            ny = impSource.getHeight();
            nz = impSource.getStackSize();

            Imax = impSource.getDisplayRangeMax();
            Imin = impSource.getDisplayRangeMin();
            if (type == ImagePlus.COLOR_RGB) {
                if (nz > 1) {
                    IJ.showMessage("Note", "This version of the plugin does not handle color stacks.");
                    return;
                } else {
                    COLOR = true;
                    nz = 3;
                    original = Builder.create(nx, ny, nz, ImageWare.DOUBLE);
                    ImageWare[] temp = Builder.createColors(impSource);
                    original.putXY(0, 0, 0, temp[0]);
                    original.putXY(0, 0, 1, temp[1]);
                    original.putXY(0, 0, 2, temp[2]);
                }
            } else {
                COLOR = false;
                original = Builder.create(impSource);
            }
            /*
            if (Math.IEEEremainder(nx, 2.0) != 0 && Math.IEEEremainder(ny, 2.0) != 0) {
                IJ.log("Your data should have at least one even dimension.");
                return;
            }
             * 
             */
            if (nx < Nmin || ny < Nmin) {
                IJ.log("The size of your data is inapropriate.");
                return;
            }
            if (original == null) {
                IJ.log("Unable to create the data set.");
                return;
            }
            IJ.showStatus("Denoising in progress...");
            original = original.convert(ImageWare.DOUBLE);
            nz = original.getSizeZ();
            nxe = (int) (Math.ceil((double) nx / Nmin) * Nmin);
            nye = (int) (Math.ceil((double) ny / Nmin) * Nmin);
            if (nxe != nx || nye != ny) {
                original = Operations.symextend2D(original, nxe, nye, Ext);
            } else {
                Ext[0] = 0;
                Ext[1] = 0;
            }
            AlphaHat = new double[nz];
            DeltaHat = new double[nz];
            SigmaHat = new double[nz];
            denoising = new Denoising(original, AlphaHat, DeltaHat, SigmaHat, FRAMEWISE, CS, NBFRAME);
            denoising.setLog(LOG);

            String argsEstimation[] = split(Macro.getValue(options, "estimation", "Auto Global"));
            for (int u = 0; u < argsEstimation.length; u++) {
                IJ.log("argsEstimation["+u+"]: "+argsEstimation[u]);
            }
            if (argsEstimation.length != 2 && argsEstimation.length != 4) {
                IJ.error("The estimation parameters is incorrect. Correct example:\"estimation=Auto Global\" or \"estimation=Manual 30.0 3.0 40.0\" ");
                return;
            }

            if (argsEstimation[0].toLowerCase().equals("auto")) {
                denoising.setFramewise(argsEstimation[1].toLowerCase().equals("individual"));
                //System.out.print(argsEstimation[1].toLowerCase().equals("individual"));
                denoising.estimateNoiseParameters();
            } else if (argsEstimation[0].toLowerCase().equals("manual")) {
                Alpha = (new Double(argsEstimation[1])).doubleValue();
                Delta = (new Double(argsEstimation[2])).doubleValue();
                Sigma = (new Double(argsEstimation[3])).doubleValue();
                for (int i = 0; i < nz; i++) {
                    AlphaHat[i] = Alpha;
                    DeltaHat[i] = Delta;
                    SigmaHat[i] = Sigma;
                }
            } else {
                IJ.error("The estimation parameters is incorrect. Correct example:\"estimation=Auto Global\" or \"estimation=Manual 30.0 3.0 40.0\" ");
                return;
            }

            String argsParameters[] = split(Macro.getValue(options, "parameters", "3 4"));
            for (int u = 0; u < argsParameters.length; u++) {
            	IJ.log("argsParameters["+u+"]: "+argsParameters[u]);
            }
            if (argsParameters.length != 2) {
                IJ.error("The parameters is incorrect. Correct example:\"parameters=3 4\" ");
                return;
            }
            NBFRAME = Math.max(1, Math.min(nz, (int) (new Double(argsParameters[0])).doubleValue()));
            CS = Math.max(1, Math.min(10, (int) (new Double(argsParameters[1])).doubleValue()));
            denoising.setCycleSpins(CS);
            denoising.setMultiFrame(NBFRAME);
            denoising.perform();
            output = denoising.getOutput();
            if (nxe != nx || nye != ny) {
                output = Operations.crop2D(output, nx, ny, Ext);
            }
            display();
            //output.show("Denoised " + impSource.getTitle());
            //----------------------------------------------------------------------------------
            IJ.log("-------------- SUMMARY --------------");
            IJ.log("Noise parameters used for denoising: \"" + titleIn + "\"");
            for (int i = 0; i < nz; i++) {
                IJ.log("Frame " + (i + 1) + ": Alpha = " + IJ.d2s(AlphaHat[i], 3) + " Delta = " + IJ.d2s(DeltaHat[i], 3) + " Sigma = " + IJ.d2s(SigmaHat[i], 3));
            }
            IJ.log("Number of adjacent frames = " + NBFRAME);
            IJ.log("Number of cycle-spin(s) = " + CS);
            IJ.log("Maximum number of concurrent threads: " + denoising.getMaxThread());
            IJ.log("The whole processing required " + IJ.d2s((System.currentTimeMillis() - time) / 1000.0, 2) + " s.");
            IJ.log("-------------------------------------");
        } else {
            new DenoiseDialog();
        }
    }

    public synchronized void display() {
        IJ.showStatus("Displaying the results...");
        if (COLOR) {
            ImageWare[] temp = new ImageWare[3];
            for (int i = 0; i < 3; i++) {
                temp[i] = Builder.create(nx, ny, 1, ImageWare.DOUBLE);
                output.getXY(0, 0, i, temp[i]);
                temp[i].rescale(0, 255);
            }
            Display.showColor(titleOut, temp[0], temp[1], temp[2]);
        } else {
            output.show(titleOut);
            ImagePlus impOutput = WindowManager.getCurrentImage();
            impOutput.setDisplayRange(Imin, Imax);
            ColorModel cm = new IndexColorModel(8, 256, red, green, blue);
            ImageProcessor ip = impOutput.getProcessor();
            ip.setColorModel(cm);
            impOutput.updateImage();
        }
        notify();
    }

    private String[] split(String s) {
        StringTokenizer t = new StringTokenizer(s);
        String[] items = new String[t.countTokens()];
        for (int k = 0; (k < items.length); k++) {
            items[k] = t.nextToken();
        }
        return (items);
    }
}