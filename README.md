# RNA_Scope_Nerve

* **Developed for:** Vianney
* **Team:** Brunet
* **Date:** March 2023
* **Software:** Fiji


### Images description

3D images taken with a x63 objective on an Airyscan

2 channels:
  1. *488:* Gene1 dots
  2. *555:* Gene2 dots
  
A *.roi* or *.zip* file containing ROI(s) can be provided with each image.

### Plugin description

In each ROI,
* Detect Gene1 dots with Median filtering + DoG filtering + MaxEntropy thresholding
* Detect Gene2 dots with Median filtering + DoG filtering + MaxEntropy thresholding
* Estimate number of genes with foci Single foci estimated volume

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin

### Version history

Version 1 released on March 2, 2023.
