=== Do NOT run this as one code. There are three sections that require other steps be taken in between. ===

1a) Compute ARGUS camera extrinsics

1b) Define timestacks and get .pix file for i2R system

------ camera acquires timestacks -----

2) extract timestacks from i2R product

------ run segmentation_gym on images -----
   
3a) extract horizontal runup line from segmented image

3b) project horizontal runup line onto MOPS vertical elevation

======= i2Rgus =====

.pix files were placed in arguseyes/build folder
and .xml file is what is executed by the crontab

======= Segmentation Gym =====

Here's the main model
Doodleverse/segmentation_gym: https://github.com/Doodleverse/segmentation_gym
 
Here's the newer segformer model:
Doodleverse/Segmentation Gym SegFormer models for 2-class (water, other) segmentation of greyscale CoastCam runup timestack imagery (zenodo.org)
https://zenodo.org/records/11167477
 
And here's the 'old' resunet model (this is the one that I can get running and doesn't require an nvidia graphics card)
Doodleverse/Segmentation Gym Res-UNet models for 2-class (water, other) segmentation of CoastCam runup timestack imagery (zenodo.org)
https://zenodo.org/records/7921971
 
I had errors for lib- something before and running this at the beginning of the process helped:
sudo apt update && \
sudo apt install libgl1-mesa-dri libegl1 libglu1-mesa libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2-data libasound2-plugins libxi6 libxtst6
 
