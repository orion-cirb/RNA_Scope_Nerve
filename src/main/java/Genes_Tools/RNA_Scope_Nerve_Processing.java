package Genes_Tools;




import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.plugin.RGBStackMerge;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import org.apache.commons.io.FilenameUtils;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;



/**
 *
 * @author phm
 */

public class RNA_Scope_Nerve_Processing {
    
    public CLIJ2 clij2 = CLIJ2.getInstance();
    
    // min size for dots
    public double minFoci = 0.05;
    public double maxFoci = 10;
    public double singleDotVol = 0.06;

    private double minDOGDots = 1;
    private double maxDOGDots = 2;
    private String geneThreshold = "MaxEntropy";
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    public Calibration cal = new Calibration();    
    public float pixVol = 0;
    
     /**
     * check  installed modules
     * @return 
     */
    public boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }

        
   /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No Image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
    
    /**
     * Dialog
     */
    public String[] dialog(String[] channels) {
        String[] channelsName = {"Gene1 : ", "Gene2 : "};
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 20, 0);
        gd.addImage(icon);
        gd.addMessage("Channels selection", Font.getFont("Monospace"), Color.blue);
        for (int n = 0; n < channelsName.length; n++) {
            gd.addChoice(channelsName[n], channels, channels[0]);
        }
        gd.addMessage("Dots filter", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min foci volume : ", minFoci, 2, 6, "µm3");
        gd.addNumericField("Max foci volume : ", maxFoci, 2, 6, "µm3");
        gd.addNumericField("Single foci estimated volume : ", singleDotVol, 2, 6, "µm3");
        // Calibration
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY pixel size : ", cal.pixelWidth, 3);
        gd.addNumericField("Z pixel size : ", cal.pixelDepth, 3);
        gd.showDialog();
        String[] ch = new String[channelsName.length];
        for (int i = 0; i < channelsName.length; i++)
            ch[i] = gd.getNextChoice();

        minFoci = gd.getNextNumber();
        maxFoci= gd.getNextNumber();
        singleDotVol = gd.getNextNumber();
        cal.pixelWidth = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();        
        pixVol = (float) (cal.pixelWidth*cal.pixelHeight*cal.pixelDepth);
        if (gd.wasCanceled())
                ch = null;
        return(ch);
    } 
    
    /**
     * Find channels name
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs+1];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelFluor(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break;    
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break; 
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        channels[chs] = "None";
        return(channels);     
    }
    
    /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public Calibration findImageCalib(IMetadata meta) {
        // read image calibration
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
        return(cal);
    }
    
    
    /**
     *
     * @param img
     */
    public void closeImages(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    
    /**
     * Difference of Gaussians 
     * Using CLIJ2
     * @param imgCL
     * @param sizeX1
     * @param sizeY1
     * @param sizeX2
     * @param sizeY2

     * @return imgGauss
     */ 
    public ClearCLBuffer DOG(ClearCLBuffer imgCL, double sizeX1, double sizeY1, double sizeZ1, double sizeX2, double sizeY2, double sizeZ2) {
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian3D(imgCL, imgCLDOG, sizeX1, sizeY1, sizeZ1, sizeX2, sizeY2, sizeZ2);
        clij2.release(imgCL);
        return(imgCLDOG);
    }
    
    
    /**
     * Threshold 
     * USING CLIJ2
     * @param imgCL
     * @param thMed
     * @param fill 
     */
    public ClearCLBuffer threshold(ClearCLBuffer imgCL, String thMed) {
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        clij2.release(imgCL);
        return(imgCLBin);
    }
    
    /**
     * return objects population in an binary image
     * @param img
     * @return pop
     */
    public Objects3DIntPopulation getPopFromImage(ImagePlus img, Calibration cal) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        labels.setCalibration(cal);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        return pop;
    }
    
    
    
     /**
     * Remove object with size < min and size > max
     * @param pop
     * @param min
     * @param max
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
    
     /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public void clearOutSide(ImagePlus img, Roi roi) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), Roi.FREEROI);
        poly.setLocation(0, 0);
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(poly);
            ip.setBackgroundValue(0);
            ip.fillOutside(poly);
        }
        img.deleteRoi();
        img.updateAndDraw();
    }

    /**
     * Find genes population
     * @param imgGene
     * @return genePop
     */
    public Objects3DIntPopulation findGenesPop(ImagePlus imgGene, Roi roi) {
        IJ.showStatus("Finding gene dots ...");
        ClearCLBuffer imgCLMed = clij2.push(imgGene);
        ClearCLBuffer imgCLDOG = DOG(imgCLMed, minDOGDots, minDOGDots, minDOGDots, maxDOGDots, maxDOGDots, maxDOGDots);
        clij2.release(imgCLMed);
        ClearCLBuffer imgCLBin = threshold(imgCLDOG, "MaxEntropy"); 
        clij2.release(imgCLDOG);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        if (roi != null)
            clearOutSide(imgBin, roi);
        Objects3DIntPopulation genePop = getPopFromImage(imgBin, cal);
        popFilterSize(genePop, minFoci, maxFoci);
        closeImages(imgBin);
        clij2.release(imgCLBin);       
        return(genePop);
    }

    /**
     * Find sum volume of objects  
     * @param dotsPop
     * @return vol
     */
    
    public double findPopVolume (Objects3DIntPopulation dotsPop) {
        IJ.showStatus("Findind object's volume");
        double sumVol = 0;
        for(Object3DInt obj : dotsPop.getObjects3DInt()) {
            sumVol += new MeasureVolume(obj).getVolumeUnit();
        }
        return(sumVol);
    }
    
    
   /**
     * Find roi volume
     * @param roi
     * @param img
     * @return volume
     */
    public double roiVolume(Roi roi, ImagePlus imgAstro) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), PolygonRoi.FREEROI); 
        poly.setLocation(0, 0);
        ImageProcessor ip = imgAstro.getProcessor();
        ip.setRoi(poly);
        ImageStatistics stats = ImageStatistics.getStatistics(ip, Measurements.AREA, cal);
        double volume = stats.area * imgAstro.getNSlices();
        return(volume);
    }
    
    
    /**
     * save images objects population
     * @param img
     * @param gene1Pop
     * @param gene2Pop
     * @param path
     */
    public void saveGenesImage (ImagePlus img, Objects3DIntPopulation gene1Pop, Objects3DIntPopulation gene2Pop, String path) {
        // red gene1 , green gene2
        ImageHandler imhGene1 = ImageHandler.wrap(img).createSameDimensions();
        ImageHandler imhGene2 = (gene2Pop.getNbObjects() > 0) ? imhGene1.duplicate() : null;
        // draw genes population
        for (Object3DInt ob : gene1Pop.getObjects3DInt())
            ob.drawObject(imhGene1, 255);
        if (imhGene2 != null)
            for (Object3DInt ob : gene2Pop.getObjects3DInt())
                ob.drawObject(imhGene2, 255);
        ImagePlus[] imgColors = {imhGene1.getImagePlus(), imhGene2.getImagePlus()};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");

        // Save images
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(path);
        imhGene1.closeImagePlus();
        if (imhGene2 != null)
            imhGene2.closeImagePlus();
    }
   
}
