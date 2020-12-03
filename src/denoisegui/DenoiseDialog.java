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

import ij.IJ;
import ij.ImagePlus;
import ij.LookUpTable;
import ij.WindowManager;
import ij.gui.GUI;
import ij.process.ImageProcessor;
import imageware.Builder;
import imageware.Display;
import imageware.ImageWare;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.TextEvent;
import java.awt.event.TextListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.awt.image.ColorModel;
import java.awt.image.IndexColorModel;
import java.util.Hashtable;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JRadioButton;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import additionaluserinterface.GridPanel;
import additionaluserinterface.WalkBar;
import denoise.Denoising;
import denoise.Operations;

public class DenoiseDialog extends JDialog
	implements ChangeListener, ActionListener, ItemListener, WindowListener, TextListener, Runnable
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String  defaultMessage  = "(c) 2010 EPFL, BIG";
    private WalkBar walk            = new WalkBar(defaultMessage, true, false, true);
    private Thread  thread          = null;

    private double 		Alpha = 1.0, Delta = 0.0, Sigma = 0.0;
    private double[]    AlphaHat, DeltaHat, SigmaHat;	// Estimated gain & offset of the detectors, AGWN std
    private double 		Imax = 0.0, Imin = 0.0;
    private int         CS = 4;
    private int         NBFRAME = 3;
    private int         Nmin = 16;
    private Denoising   denoising;

    private int nbMaxCharInput = 25;

    private ImagePlus impSource  = null;
    private ImagePlus impOutput  = null;

    private String titleInput  = "";

    private ImageWare original = null;		// Input image
    private ImageWare output = null;		// Output image
    private boolean   COLOR = false;	    // Index for color images identification
    private boolean   NOISEPARAMS = false;


    private int 	nx = 0;
    private int 	ny = 0;
    private int 	nz = 0;
    private int 	nxe = 0;
    private int 	nye = 0;
    private int[] 	Ext = new int[2];
    private byte[] 	red, green, blue;

    // Dialog box components
    private JButton 		bnRun = new JButton("Start");
    private JCheckBox       checkLog = new JCheckBox("Display Log");
    private JRadioButton	checkAutomatic;
    private JRadioButton 	checkManual;
    private JComboBox   	cmbAutomatic = new JComboBox(new String[] {"Global", "Individual"});
    private ButtonGroup		cbgAutomatic;
    private JTextField 		txtOutput;
    private JTextField		txtSigma;
    private JTextField		txtAlpha;
    private JTextField		txtDelta;
    private JTextField		txtnbFrame = new JTextField("",6);
    private JTextField		txtCS = new JTextField("",6);
    private JLabel 			lblCS = new JLabel("<html>Cycle-spins</html>");
    private JLabel			lblnbFrame = new JLabel("<html>Multiframe</html>");
    private JTextField  	txtTitle;
    private JLabel 			lblSigma;
    private JLabel 			lblAlpha;
    private JLabel 			lblDelta;
    private JSlider			sldCyclespin = new JSlider(JSlider.HORIZONTAL, 1, 10, CS);
    private JSlider			sldMultiframe = new JSlider(JSlider.HORIZONTAL, 1, 11, NBFRAME);
    private Timer       	timer;
	
    /**
    * Constructor creates the dialog box.
    */

    public DenoiseDialog() {
        super(new Frame(), "PureDenoise");
        walk.fillAbout(
                "PureDenoise",
                "Version 12/06/2010",
                "Biomedical Image Denoising",
                "Florian Luisier",
                "Biomedical Imaging Group (BIG)<br>Ecole Polytechnique F&eacute;d&eacute;rale de Lausanne (EPFL)<br>Lausanne, Switzerland",
                "2010",
                "<p style=\"text-align:left\"><b>Full Info:</b>"+
                "<br>http://bigwww.epfl.ch/algorithms/denoise/<br>"+
                "<br><b>References:</b>"+
                //"<small>" +
                "<br>[1] F. Luisier, C. Vonesch, T. Blu, M. Unser, Fast Interscale Wavelet Denoising of Poisson-corrupted Images, Signal Processing, vol. 90, no. 2, pp. 415-427, February 2010." +
                "<br>[2] F. Luisier, C. Vonesch, T. Blu, M. Unser, Fast Haar-Wavelet Denoising of Multidimensional Fluorescence Microscopy Data, Proceedings of the Sixth IEEE International Symposium on Biomedical Imaging: From Nano to Macro (ISBI'09), Boston MA, USA, June 28-July 1, 2009, pp. 310-313." +
                "<br>[3] F. Luisier, The SURE-LET Approach to Image Denoising, Swiss Federal Institute of Technology Lausanne, EPFL Thesis no. 4566 (2010), 232 p., January 8, 2010.<br>"+
                //"</small>" +
                "<br><b>Acknowledgements:</b>"+
                //"<small>" +
                "<br>Prof. Thierry Blu" +
                "<br>Prof. Michael Unser" +
                "<br>Dr. Daniel Sage"+
                "<br>Dr. C&eacute;dric Vonesch"
                //"</small>"
                );

        Hashtable<Integer, JLabel> labelTable1 = new Hashtable<Integer, JLabel>();
        labelTable1.put(new Integer(1), new JLabel("Fast") );
        labelTable1.put(new Integer(9), new JLabel("HQ") );

        sldCyclespin.setPaintLabels(true);
        sldCyclespin.setPaintTicks(true);
        sldCyclespin.setSnapToTicks(true);
        sldCyclespin.setLabelTable(labelTable1);
        sldCyclespin.setMinimumSize(new java.awt.Dimension(150,50));
        sldCyclespin.setPreferredSize(new java.awt.Dimension(150,50));
        sldCyclespin.setMinorTickSpacing(1);

        sldMultiframe.setPaintLabels(true);
        sldMultiframe.setSnapToTicks(true);
        sldMultiframe.setMinorTickSpacing(2);
        sldMultiframe.setPreferredSize(new Dimension(150,50));
        sldMultiframe.setMinimumSize(new java.awt.Dimension(150,50));
        sldMultiframe.setPaintTicks(true);
        sldMultiframe.setSnapToTicks(true);

        // Text Title
        txtTitle = new JTextField("", nbMaxCharInput);
        txtTitle.setEditable(false);
        txtTitle.setEnabled(false);

        // Panel Control
        GridPanel pnControls = new GridPanel();
        int row = 0;
        txtOutput = new JTextField("", nbMaxCharInput);
        txtOutput.setBackground(Color.white);
        pnControls.place(0, 0, new JLabel("Input"));
        pnControls.place(0, 1, txtTitle);
        pnControls.place(1, 0, new JLabel("Output"));
        pnControls.place(1, 1, txtOutput);

        cbgAutomatic = new ButtonGroup();
        checkAutomatic = new JRadioButton("Automatic", true);
        checkManual = new JRadioButton("Manual", false);
        cbgAutomatic.add(checkAutomatic);
        cbgAutomatic.add(checkManual);

        // Detector Gain Alpha
        txtAlpha = new JTextField("1.00  ", 5);
        txtAlpha.setForeground(Color.red);
        lblAlpha = new JLabel("Detector gain");
        lblAlpha.setForeground(Color.black);

        // Detector Offset Delta
        txtDelta = new JTextField("0.00  ", 5);
        txtDelta.setForeground(Color.red);
        lblDelta = new JLabel("Detector offset");
        lblDelta.setForeground(Color.black);

        // AGWN Std
        txtSigma = new JTextField("0.00 ", 5);
        lblSigma = new JLabel("<html>Standard deviation<br>of Gaussian noise</html>");
        txtSigma.setForeground(Color.red);
        lblSigma.setForeground(Color.black);
		
        // Panel Automatic Parameters Estimation
        GridPanel pnAutomatic = new GridPanel(true);
        pnAutomatic.place(0, 1, cmbAutomatic);
		
        // Panel for Noise Model Parameters
        GridPanel pnValue = new GridPanel(true);
        pnValue.place(0, 0, lblAlpha);
        pnValue.place(0, 1, txtAlpha);
        pnValue.place(1, 0, lblDelta);
        pnValue.place(1, 1, txtDelta);
        pnValue.place(2, 0, lblSigma);
        pnValue.place(2, 1, txtSigma);
		
		
        // Panel Automatic Parameters Estimation
        GridPanel pnNoise = new GridPanel("Noise estimation");
        pnNoise.place(1, 0, checkAutomatic);
        pnNoise.place(1, 1, pnAutomatic);
        pnNoise.place(3, 0, checkManual);
        pnNoise.place(3, 1, pnValue);
      

        // Panel for remaining parameters
        txtCS.setForeground(Color.red);
        txtnbFrame.setForeground(Color.red);
        txtCS.setText("4 cycles");
        txtnbFrame.setText("3 frames");
        txtCS.setEditable(false);
        txtnbFrame.setEditable(false);
        GridPanel pnOther = new GridPanel("Denoising parameters");
        pnOther.place(0, 0, lblCS);
        pnOther.place(0, 1, sldCyclespin);
        pnOther.place(0, 2, txtCS);
        pnOther.place(1, 0, lblnbFrame);
        pnOther.place(1, 1, sldMultiframe);
        pnOther.place(1, 2, txtnbFrame);

        // Panel Buttons
        GridPanel pnButtons = new GridPanel(false);
        pnButtons.place(0, 0, 1, 1, 3, checkLog);
        pnButtons.place(0, 1, 1, 1, 3, bnRun);

        // Panel Main
        GridPanel pnMain = new GridPanel(false, 7);
        row = 0;
        pnMain.place(row++, 0, pnControls);
        pnMain.place(row++, 0, pnNoise);
        pnMain.place(row++, 0, pnOther);
        pnMain.place(row++, 0, pnButtons);
        pnMain.place(row++, 0, walk);

        // Add Listeners
        walk.getButtonClose().addActionListener(this);
        bnRun.addActionListener(this);
        checkAutomatic.addItemListener(this);
        checkManual.addItemListener(this);
        cmbAutomatic.addItemListener(this);
        checkLog.addItemListener(this);
        sldCyclespin.addChangeListener(this);
        sldMultiframe.addChangeListener(this);
        addWindowListener(this);

        bnRun.setEnabled(false);
        checkAutomatic.setEnabled(false);
        checkManual.setEnabled(false);
        cmbAutomatic.setEnabled(false);
        checkLog.setEnabled(false);
        txtCS.setEditable(false);
        txtCS.setForeground(Color.gray);
        txtnbFrame.setEditable(false);
        txtnbFrame.setForeground(Color.gray);
        txtAlpha.setEditable(false);
        txtAlpha.setForeground(Color.gray);
        lblAlpha.setForeground(Color.gray);
        txtSigma.setEditable(false);
        txtSigma.setForeground(Color.gray);
        lblSigma.setForeground(Color.gray);
        txtDelta.setEditable(false);
        txtDelta.setForeground(Color.gray);
        lblDelta.setForeground(Color.gray);

        sldCyclespin.setEnabled(false);
        sldMultiframe.setEnabled(false);
        
        // Building the main panel
        add(pnMain);
        pack();
        GUI.center(this);
        setVisible(true);
        IJ.wait(250); 	// work around for Sun/WinNT bug

    }

    /**
    * Prepare the image sequence for further use, deals with gray-scale stacks.
    *
    * The image type is identified and the appropriate ImageWare object is constructed.
    */
    private void selectInputImage()
    {
        bnRun.setEnabled(false);
        checkAutomatic.setEnabled(false);
        checkManual.setEnabled(false);
        cmbAutomatic.setEnabled(false);
        checkLog.setEnabled(false);
        txtCS.setEditable(false);
        txtCS.setForeground(Color.gray);
        txtnbFrame.setEditable(false);
        txtnbFrame.setForeground(Color.gray);
        txtAlpha.setEditable(false);
        txtAlpha.setForeground(Color.gray);
        lblAlpha.setForeground(Color.gray);
        txtSigma.setEditable(false);
        txtSigma.setForeground(Color.gray);
        lblSigma.setForeground(Color.gray);
        txtDelta.setEditable(false);
        txtDelta.setForeground(Color.gray);
        lblDelta.setForeground(Color.gray);
  
        txtOutput.setEditable(false);
        txtOutput.setEnabled(false);
        sldCyclespin.setEnabled(false);
        sldMultiframe.setEnabled(false);

        if(thread!=null){
            bnRun.setEnabled(true);
            return;
        }
        impSource = WindowManager.getCurrentImage();
        if (impSource == null) {
            txtTitle.setText("Please open an image.");
            txtOutput.setText("Please open an image.");
            thread = null;
            return;
        }
        LookUpTable lut = impSource.createLut();
        red   = lut.getReds();
        green = lut.getGreens();
        blue  = lut.getBlues();

        nx = impSource.getWidth();
        ny = impSource.getHeight();
        nz = impSource.getStackSize();

        if (nx<Nmin || ny<Nmin) {
            txtTitle.setText("The size of your data is inapropriate.");
            txtOutput.setText("The size of your data is inapropriate.");
            thread = null;
            return;
        }

        Imax = impSource.getDisplayRangeMax();
        Imin = impSource.getDisplayRangeMin();

        int type = impSource.getType();
        if(impSource.isComposite()){
            txtTitle.setText("This plugin does not handle composite.");
            txtOutput.setText("This plugin does not handle composite.");
            thread = null;
            return;
        }
        if(type==ImagePlus.COLOR_RGB) {
            if(nz>1){
                txtTitle.setText("This plugin does not handle color stacks.");
                txtOutput.setText("This plugin does not handle color stacks.");
                thread = null;
                return;
            }
            else{
                COLOR    = true;
                nz       = 3;
                original = Builder.create(nx,ny,nz,ImageWare.DOUBLE);
                ImageWare[] temp = Builder.createColors(impSource);
                original.putXY(0,0,0,temp[0]);
                original.putXY(0,0,1,temp[1]);
                original.putXY(0,0,2,temp[2]);
            }

        }
        else {
            COLOR    = false;
            original = Builder.create(impSource);
        }
        if (original == null) {
            txtTitle.setText("Unable to create the data set.");
            txtOutput.setText("Unable to create the data set.");
            thread = null;
            return;
        }
        bnRun.setText("Start");
        original = original.convert(ImageWare.DOUBLE);
        nz       = original.getSizeZ();     
        nxe = (int)(Math.ceil((double)nx/Nmin)*Nmin);
        nye = (int)(Math.ceil((double)ny/Nmin)*Nmin);
        if(nxe!=nx||nye!=ny)
            original = Operations.symextend2D(original,nxe,nye,Ext);
        else{
            Ext[0] = 0;
            Ext[1] = 0;
        }
        
        if(AlphaHat==null || AlphaHat.length!=nz || !impSource.getTitle().matches(titleInput)){
            AlphaHat 	= new double[nz];
            DeltaHat 	= new double[nz];
            SigmaHat 	= new double[nz];
            NOISEPARAMS = false;
        }
        denoising  = new Denoising(original,AlphaHat,DeltaHat,SigmaHat,checkAutomatic.isSelected(),CS,NBFRAME);
        titleInput = impSource.getTitle();
        String title = titleInput;
        if (title.length() > nbMaxCharInput-1) {
            int len = title.length();
            title = title.substring(0,nbMaxCharInput-6) + "..." + title.substring(len-3, len);
        }
        sldCyclespin.setEnabled(true);
        if (nz == 1) {
            sldMultiframe.setEnabled(false);
            sldMultiframe.setMaximum(1);
        }
        else {
            sldMultiframe.setEnabled(true);
            sldMultiframe.setMaximum(Math.min(11, nz));
            sldMultiframe.setValue(NBFRAME);
        }

        txtTitle.setText(title);
        txtTitle.setCaretPosition(0);
        bnRun.setEnabled(true);
        txtOutput.setText("Denoised-"+impSource.getTitle());
        txtOutput.setCaretPosition(0);
        txtOutput.setEditable(true);
        txtOutput.setEnabled(true);

        Alpha = getDoubleFromTextField(txtAlpha);
        Delta = getDoubleFromTextField(txtDelta);
        Sigma = getDoubleFromTextField(txtSigma);
        if(checkAutomatic.isSelected() && (cmbAutomatic.getSelectedIndex()==0) &&
           (Math.abs(Alpha-AlphaHat[0])>1e-2||
            Math.abs(Sigma-SigmaHat[0])>1e-2||
            Math.abs(Delta-DeltaHat[0])>1e-2)){
            denoising.setFramewise(false);
            NOISEPARAMS = denoising.estimateNoiseParameters();
            AlphaHat = denoising.getAlpha();
            DeltaHat = denoising.getDelta();
            SigmaHat = denoising.getSigma();
            Alpha = AlphaHat[0];
            Delta = DeltaHat[0];
            Sigma = SigmaHat[0];
            txtAlpha.setText(IJ.d2s(Alpha,2));
            txtDelta.setText(IJ.d2s(Delta,2));
            txtSigma.setText(IJ.d2s(Sigma,2));
        }

        if(checkAutomatic.isSelected()){
            cmbAutomatic.setEnabled(true);
            denoising.setFramewise(cmbAutomatic.getSelectedIndex()==1);
        }
        else{
            txtAlpha.setEditable(true);
            txtAlpha.setForeground(Color.red);
            lblAlpha.setForeground(Color.black);
            txtSigma.setEditable(true);

            txtSigma.setForeground(Color.red);
            lblSigma.setForeground(Color.black);
            txtDelta.setEditable(true);

            txtDelta.setForeground(Color.red);
            lblDelta.setForeground(Color.black);
            txtnbFrame.setEditable(true);
            txtnbFrame.setForeground(Color.red);
            denoising.setFramewise(false);
            
            NOISEPARAMS = true;
        }
        bnRun.setEnabled(true);
        checkAutomatic.setEnabled(true);
        checkManual.setEnabled(true);
        txtCS.setForeground(Color.red);
        txtnbFrame.setForeground(Color.red);
        checkLog.setEnabled(true);
    }


    /**
    * Implements the textValueChanged for the TexttListener.
    */
    public synchronized void textValueChanged(TextEvent e)
    {
        notify();
    }

    /**
    * Implements the adjustmentValueChanged for the ChangeListener.
    */
    public synchronized void stateChanged(ChangeEvent e)
    {
        if (e.getSource() == sldCyclespin) {
            CS = (int)sldCyclespin.getValue();
            if (sldCyclespin.getValue() == 1)
                txtCS.setText("1 cycle  ");
            else
                txtCS.setText("" + sldCyclespin.getValue() + " cycles");
        }
        if (e.getSource() == sldMultiframe) {
            NBFRAME = (int)sldMultiframe.getValue();
            if (sldMultiframe.getValue() == 1)
                txtnbFrame.setText("1 frame   ");
            else
                txtnbFrame.setText("" + sldMultiframe.getValue() + " frames");
        }
        notify();
    }

    /**
    *  Actualize the dialog box in case of Automatic checkbox state changes.
    */
    public synchronized void itemStateChanged(ItemEvent e)
    {
        if(checkAutomatic.isSelected()){
            cmbAutomatic.setEnabled(true);
            txtAlpha.setEditable(false);
            txtAlpha.setForeground(Color.gray);
            lblAlpha.setForeground(Color.gray);
            txtSigma.setEditable(false);
            txtSigma.setForeground(Color.gray);
            lblSigma.setForeground(Color.gray);
            txtDelta.setEditable(false);
            txtDelta.setForeground(Color.gray);
            lblDelta.setForeground(Color.gray);
            denoising.setFramewise(cmbAutomatic.getSelectedIndex()==1);
            if (cmbAutomatic.getSelectedIndex()==0 &&
                (Math.abs(Alpha - AlphaHat[0]) > 1e-2
              || Math.abs(Sigma - SigmaHat[0]) > 1e-2
              || Math.abs(Delta - DeltaHat[0]) > 1e-2)) {
                checkAutomatic.setEnabled(false);
                checkManual.setEnabled(false);
                denoising.setFramewise(false);
                NOISEPARAMS = denoising.estimateNoiseParameters();
                AlphaHat = denoising.getAlpha();
                DeltaHat = denoising.getDelta();
                SigmaHat = denoising.getSigma();
                txtAlpha.setText(IJ.d2s(AlphaHat[0], 2));
                txtDelta.setText(IJ.d2s(DeltaHat[0], 2));
                txtSigma.setText(IJ.d2s(SigmaHat[0], 2));
            }
        }
        else{
            cmbAutomatic.setEnabled(false);
            txtAlpha.setEditable(true);
            txtAlpha.setForeground(Color.red);
            lblAlpha.setForeground(Color.black);
            txtSigma.setEditable(true);
            txtSigma.setForeground(Color.red);
            lblSigma.setForeground(Color.black);
            txtDelta.setEditable(true);
            txtDelta.setForeground(Color.red);
            lblDelta.setForeground(Color.black);
            txtnbFrame.setEditable(true);
            txtnbFrame.setForeground(Color.red);
            lblnbFrame.setForeground(Color.black);
            denoising.setFramewise(false);
            NOISEPARAMS = true;
        }
        Alpha = getDoubleFromTextField(txtAlpha);
        Delta = getDoubleFromTextField(txtDelta);
        Sigma = getDoubleFromTextField(txtSigma);
        checkAutomatic.setEnabled(true);
        checkManual.setEnabled(true);
        notify();
    }

    /**
    */
    public void run() {
        if(!NOISEPARAMS){
            IJ.error("The automatic noise parameters estimation failed!\nSet the noise parameters manually.");
            checkAutomatic.setEnabled(false);
            cmbAutomatic.setEnabled(false);
            checkManual.setEnabled(true);
            thread = null;
            return;
        }
        double time = System.currentTimeMillis();
        walk.reset();
        walk.setMessage("Denoising in progress...");
        IJ.showStatus("Initializing...");
        txtOutput.setEnabled(false);
        bnRun.setText("Stop");
        checkAutomatic.setEnabled(false);
        checkManual.setEnabled(false);
        cmbAutomatic.setEnabled(false);
        sldCyclespin.setEnabled(false);
        sldMultiframe.setEnabled(false);
        checkLog.setEnabled(false);
        txtAlpha.setEditable(false);
        txtAlpha.setForeground(Color.gray);
        lblAlpha.setForeground(Color.gray);
        txtSigma.setEditable(false);
        txtSigma.setForeground(Color.gray);
        lblSigma.setForeground(Color.gray);
        txtDelta.setEditable(false);
        txtDelta.setForeground(Color.gray);
        lblDelta.setForeground(Color.gray);
        txtCS.setEditable(false);
        txtCS.setForeground(Color.gray);
        txtnbFrame.setEditable(false);
        txtnbFrame.setForeground(Color.gray);
        boolean OK = getParameters();
        if(!OK){
            output = null;
            thread = null;
            IJ.showStatus("Process aborted.");
            walk.reset();
            walk.setMessage(defaultMessage);
            return;
        }
        setCursor(new Cursor(Cursor.WAIT_CURSOR));
        timer = new Timer(walk, denoising);
        timer.start();
        denoising.perform();
        timer.stop();
        output = denoising.getOutput();
        bnRun.setText("Start");		
        if(output==null) {
            thread = null;
            IJ.showStatus("Process aborted.");
            setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
            walk.reset();
            walk.setMessage(defaultMessage);
            return;
        }
        //----------------------------------------------------------------------------------
        IJ.showStatus("Displaying the results...");
        if (nxe != nx || nye != ny) {
            output = Operations.crop2D(output,nx,ny,Ext);
        }
        //----------------------------------------------------------------------------------
        if(COLOR){
            ImageWare[] temp = new ImageWare[3];
            for(int i=0;i<3;i++) {
                temp[i] = Builder.create(nx,ny,1,ImageWare.DOUBLE);
                output.getXY(0,0,i,temp[i]);
                temp[i].rescale(0,255);
            }
            Display.showColor(txtOutput.getText(),temp[0],temp[1],temp[2]);
        }
        else {
            output.show(txtOutput.getText());
            impOutput = WindowManager.getCurrentImage();
            impOutput.setDisplayRange(Imin,Imax);
            ColorModel cm = new IndexColorModel(8, 256, red, green, blue);
            ImageProcessor ip = impOutput.getProcessor();
            ip.setColorModel(cm);
            impOutput.updateImage();
        }
        thread = null;
        if(checkLog.isSelected()){
            IJ.log("-------------- SUMMARY --------------");
            IJ.log("Parameters used for denoising: \""+titleInput+"\"");
            for(int i=0;i<nz;i++)
                IJ.log("Frame "+(i+1)+": Alpha = "+IJ.d2s(AlphaHat[i],3)+" Delta = "+IJ.d2s(DeltaHat[i],3)+" Sigma = "+IJ.d2s(SigmaHat[i],3));
            IJ.log("Number of adjacent frames: "+NBFRAME);
            IJ.log("Number of cycle-spin(s): "+CS);
            IJ.log("Maximum number of concurrent threads: "+denoising.getMaxThread());
            IJ.log("The whole processing required "+IJ.d2s((System.currentTimeMillis()-time)/1000.0,2)+" s.");
        }
        setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
        walk.finish("End ...");
    }

    /**
    * Implements the actionPerformed for the ActionListener.
    */
    public synchronized  void actionPerformed(ActionEvent e)
    {
        if (e.getSource() == walk.getButtonClose()) {
            dispose();
        }
        else if(e.getSource()==bnRun) {
            if(thread==null) {
                thread = new Thread(this);
                thread.setPriority(Thread.MIN_PRIORITY);
                thread.start();
            }
            else{
                denoising.setStop(true);
                IJ.showStatus("Process aborted.");
                setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
                timer.stop();
                walk.reset();
                walk.setMessage(defaultMessage);
                thread = null;
                output = null;
                selectInputImage();
                return;
            }

        }
        notify();
    }

    /**
    * Implements the windowActivated method for the WindowListener.
    */
    public void windowActivated(WindowEvent e)
    {
        selectInputImage();
    }

    /**
    * Implements the windowClosing method for the WindowListener.
    */
    public void windowClosing(WindowEvent e)
    {
        dispose();
    }

    /**
    * Implements the windowClosed method for the WindowListener.
    */
    public void windowClosed(WindowEvent e)
    {
    }

    /**
    * Implements the windowDeactivated method for the WindowListener.
    */
    public void windowDeactivated(WindowEvent e)
    {
    }

    /**
    * Implements the windowDeiconified method for the WindowListener.
    */
    public void windowDeiconified(WindowEvent e)
    {
    }

    /**
    * Implements the windowIconified method for the WindowListener.
    */
    public void windowIconified(WindowEvent e)
    {
    }

    /**
    * Implements the windowOpened method for the WindowListener.
    */
    public void windowOpened(WindowEvent e)
    {
    }


    /**
    * Convert a String to a double.
    */
    final private double getDoubleFromTextField(final JTextField txt)
    {
        double d = 0.0;
        try {
            d = (new Double(txt.getText())).doubleValue();
        }

        catch (Exception e) {
            if (e instanceof NumberFormatException)
                IJ.error("Error in number format");
        }
        return d;
    }

    /**
    * Get all parameters of the dialog box.
    */
    final public boolean getParameters()
    {
        Alpha = getDoubleFromTextField(txtAlpha);
        Delta = getDoubleFromTextField(txtDelta);
        Sigma = getDoubleFromTextField(txtSigma);
        boolean framewise = cmbAutomatic.getSelectedIndex()==1;
        denoising.setFramewise(framewise);
        denoising.setLog(checkLog.isSelected());
        if(checkAutomatic.isSelected()){
            boolean OK = denoising.estimateNoiseParameters();
            if(!OK){
                if(!denoising.getStop()){
                    IJ.error("The automatic noise parameters estimation failed!\nSet the noise parameters manually.");
                    checkAutomatic.setEnabled(false);
                    cmbAutomatic.setEnabled(false);
                    checkManual.setEnabled(true);
                }
                return false;
            }
            AlphaHat = denoising.getAlpha();
            DeltaHat = denoising.getDelta();
            SigmaHat = denoising.getSigma();
            Alpha    = AlphaHat[0];
            Delta    = DeltaHat[0];
            Sigma    = SigmaHat[0];
            txtAlpha.setText(IJ.d2s(Alpha,2));
            txtDelta.setText(IJ.d2s(Delta,2));
            txtSigma.setText(IJ.d2s(Sigma,2));
        }
        if(checkManual.isSelected()){
            if(Alpha<=0){
                IJ.error("The detector gain should be strictly positive!");
                NOISEPARAMS = false;
                return false;
            }
            if(Sigma<0){
                IJ.error("The standard deviation of the AWGN should be non-negative!");
                NOISEPARAMS = false;
                return false;
            }
            for(int i=0;i<nz;i++){
                AlphaHat[i] = Alpha;
                DeltaHat[i] = Delta;
                SigmaHat[i] = Sigma;
            }
            denoising.setAlpha(AlphaHat);
            denoising.setDelta(DeltaHat);
            denoising.setSigma(SigmaHat);
        }
        NBFRAME = (int)sldMultiframe.getValue();
        CS      = (int)sldCyclespin.getValue();
        if(NBFRAME>=nz)
            NBFRAME = nz;
        else
            if(Math.IEEEremainder(NBFRAME,2.0)==0)
                NBFRAME = Math.max(1,NBFRAME-1);
        denoising.setCycleSpins(CS);
        denoising.setMultiFrame(NBFRAME);
        String txt = IJ.d2s(NBFRAME,0);
        txt = NBFRAME>1?txt+" frames":txt+" frame";
        txtnbFrame.setText(txt);

        return true;
    }

} // End of class
