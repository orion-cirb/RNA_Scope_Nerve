 RNA_Scope_Nerve

* **Developed for:** Vianney
* **Team:** Brunet
* **Date:** March 2023
* **Software:** Fiji


### Images description

3D images taken with a x63 objective on an Airyscan

3 channels:
  1. *488:* Gene1
  2. *555:* Gene2
  
A *.roi*, *.zip* or *no roi* file containing ROI(s) must be provided with each image.

### Plugin description

In each ROI:
* Detect genes 1 in 488 with DOG
* Detect genes 2 in 555 with DOG
* Compute roi volume
* Estimate number of genes with foci Single foci estimated volume

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin

### Version history

Version 1 released on March 2, 2023.
