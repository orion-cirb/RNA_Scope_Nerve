/*
 * Find dots in channel 0 and 3
 * Compute area ratio in roi
 * Author Philippe Mailly
 */



import Genes_Tools.RNA_Scope_Nerve_Processing;
import ij.*;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.awt.Rectangle;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.ArrayList;
import loci.common.Region;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;


public class RNA_Scope_Nerve implements PlugIn {

    private final boolean canceled = false;
    public String outDirResults = "";
    private String imageDir = "";

    private Genes_Tools.RNA_Scope_Nerve_Processing genes = new RNA_Scope_Nerve_Processing();
      

    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir += IJ.getDirectory("Choose images directory")+File.separator;
            if (imageDir == null) {
                return;
            }
            File inDir = new File(imageDir);
            ArrayList<String> imageFiles = genes.findImages(imageDir, "czi");
            if (imageFiles == null) {
                return;
            }
            
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find chanels, image calibration
            reader.setId(imageFiles.get(0));
            String[] channels = genes.findChannels(imageFiles.get(0), meta, reader);
            genes.cal = genes.findImageCalib(meta);
            String[] chs = genes.dialog(channels);
            if(channels == null)
                return;
            

            // create output folder
            outDirResults = inDir + File.separator+ "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            // Write headers results for results files
            FileWriter fileResults = new FileWriter(outDirResults + "results.xls", false);
            BufferedWriter outPutResults = new BufferedWriter(fileResults);
            outPutResults.write("ImageName\tROI name\tROI Volume\tGene1 Volume\tGene1 ratio volume %\tEstimated gene1 dots number\t"
                    + "Gene2 volume\tGene2 ratio volume %\tEstimated gene2 dots number\tRatio Gene1/Gene2 Nb %\tRatio Gene1/Gene2 Volume %\n");
            outPutResults.flush();            
            
            
            // Read images
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                
                // Find ROI file
                boolean roiFound = true;
                String roiFile = imageDir+rootName+".zip";
                if (!new File(roiFile).exists()) {
                    roiFile = imageDir+rootName+".roi";
                    if (!new File(roiFile).exists()) {
                        roiFound = false;
                        System.out.println("No ROI file found !");
                    }
                }
                RoiManager rm = new RoiManager(false);
                if (roiFound)
                    rm.runCommand("Open", roiFile);
                else
                    rm.add(new Roi(0, 0, reader.getSizeX(), reader.getSizeY()), 0);
                // for each roi open image and crop
                for (Roi roi : rm.getRoisAsArray()) {
                    ImporterOptions options = new ImporterOptions();
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setId(f);
                    options.setSplitChannels(true);
                    options.setCrop(true);
                    Region reg = new Region(roi.getBounds().x, roi.getBounds().y, roi.getBounds().width, roi.getBounds().height);
                    options.setCropRegion(0, reg);
                    options.doCrop();
                
                    String roiName = roi.getName();
                    Rectangle rectRoi = roi.getBounds();
                    Region regRoi = new Region(rectRoi.x, rectRoi.y, rectRoi.width, rectRoi.height);
                    options.setCropRegion(0, regRoi);
                    
                    // Open gene1
                    int indexCh = ArrayUtils.indexOf(channels, chs[0]);
                    System.out.println("Opening gene1 channel = "+ chs[0]);
                    ImagePlus imgGene1 = BF.openImagePlus(options)[indexCh];
                    double roiVol = genes.roiVolume(roi, imgGene1);
                    System.out.println("Roi "+roiName+ " vol = "+roiVol);
                    Objects3DIntPopulation gene1Pop = genes.findGenesPop(imgGene1, roi);
                    System.out.println("gene1 found = "+gene1Pop.getNbObjects());
                    double gene1Vol = genes.findPopVolume(gene1Pop);
                    double estimatedGene1Nb = (int) Math.round(gene1Vol/genes.singleDotVol);
                    
                    // Open gene2
                    ImagePlus imgGene2 = null;
                    double estimatedGene2Nb = 0;
                    double gene2Vol = 0;
                    Objects3DIntPopulation gene2Pop = new Objects3DIntPopulation();
                    if (!chs[1].equals("None")) {
                        indexCh = ArrayUtils.indexOf(channels, chs[1]);
                        System.out.println("Opening gene2 channel = "+ chs[1]);
                        imgGene2 = BF.openImagePlus(options)[indexCh];
                        // find gene2 dots
                        gene2Pop = genes.findGenesPop(imgGene2, roi);
                        System.out.println("gene2 found = "+gene2Pop.getNbObjects());
                        gene2Vol = genes.findPopVolume(gene2Pop);
                        estimatedGene2Nb = (int) Math.round(gene2Vol/genes.singleDotVol);
                        imgGene2.close();
                    }
                    
                    // Write parameters
                    IJ.showStatus("Writing parameters ...");
                    double ratioNb = (estimatedGene2Nb == 0) ?  0 : estimatedGene1Nb/estimatedGene2Nb;
                    double ratioVol = (gene2Vol == 0) ? 0 : gene1Vol/gene2Vol;
                   
                    outPutResults.write(rootName+"\t"+roiName+"\t"+roiVol+"\t"+gene1Vol+"\t"+(gene1Vol/roiVol)*100+"\t"+estimatedGene1Nb+
                            "\t"+gene2Vol+"\t"+(gene2Vol/roiVol)*100+"\t"+estimatedGene2Nb+"\t"+ratioNb*100+"\t"+ratioVol*100+"\n");
                    outPutResults.flush();
                    
                    // save image objects
                    IJ.showStatus("Save images objects ...");
                    String path = outDirResults + rootName+"_"+roiName+"_Objects.tif";
                    genes.saveGenesImage(imgGene1, gene1Pop, gene2Pop, path);  
                    imgGene1.close();
                }
            }
            outPutResults.close();
            IJ.showStatus("Process done");
        } catch (DependencyException | ServiceException | FormatException | IOException ex) {
            Logger.getLogger(RNA_Scope_Nerve.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
