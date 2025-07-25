# **CCNV**
The R package CCNV (**C**umulative **C**opy **N**umber **V**ariation) offers multiple methods for the cumulative copy number variation analysis from DNA methylation data. Two methods are integrated into the package for the processing of the data. For each method two cumulative plots, one displaying the intensity and one displaying the frequency of each aberration within the provided cohort is generated. 
### **Sample wise (SW) method**
In the sample wise method, each sample is processed consecutively after another and while this method takes longer, the user gets detailed plots giving an overview over the whole cohort and the individualities of each sample.
### **Combined Segmentation (CS) method**
In the combined segmentation method, all provided samples are processed simultaneously. This method yields quick results, giving an overview over the similarities within the whole cohort.

## **Installation Guide**
To install the CCNV package, you first have to install devtools. Therefore, open R and run the following Code snippet:
```
install.packages("devtools")
```
You can then either directly install the package from the GitHub directly using the following code snippet:
```
devtools::install_github("Neumann-Neurooncology/CCNV")
```
Or you can download the repository, navigate to the folder within R and run the installation via the following line of code:
```
devtools::install()
```

## **Usage**

After installation, multiple functions and options are offered to the user. For further details please refer to the vignette. An example of the Input table for the *cumul.CNV* function is provided below. Please note, that the naming of the *Basename*- and the *ArrayType*- column needs to be exact.

| Basename       	| ArrayType 	| OptionalAnnotation_1 	| ... 	|
|----------------	|-----------	|----------------------	|-----	|
| Path/To/Idat_1 	| EPIC      	| ...                  	| ... 	|
| Path/To/Idat_2 	| EPIC      	| ...                  	| ... 	|
| ...            	| ...       	| ...                  	| ... 	|

## Examples
Generating cumulative _Frequency_ and _Intensity Plots_, utilising our enhanced combined segmenatation approach using only one line of code.
```r
#read in sample sheet
df_inputIdat <- read.csv("path/to/file")

#call function
CCNV::cumul.CNV(df_inputIdat, segmentationMode = "multi")
```

Generating a cumulative _Frequency Plot_ using four other packages and at minimum 13 lines of code. (_Intensity Plots_ generally not supported.)
```r
#read in sample sheet
df_inputIdat <- read.csv("path/to/file")

#reading and processing data as described using minfi
target_rgset <- minfi::read.metharray.exp(base = NULL, targets = df_inputIdat, force = TRUE)
mSetSq_f <- minfi::preprocessFunnorm(target_rgset)
target_mset <- minfi::preprocessIllumina(target_rgset)
target_mset_mapped <- minfi::mapToGenome(target_mset)

#loading data to conumee2
target_mset_loaded <- conumee2::CNV.load(target_mset) 

#generate annotation
anno <- conumee2::CNV.create_anno(array_type = "450k")  

# take public control samples
control_mset <- minfiData::MsetEx
control_mset_loaded <- conumee::CNV.load(control_mset)

# find overlapping probes between arraydata and annotations
anno@probes <-
  IRanges::subsetByOverlaps(anno@probes, granges(target_mset_mapped))


#run segmentation
x <- conumee2::CNV.fit(query =target_mset_loaded, ref = control_mset_loaded, anno)
x <- conumee2::CNV.bin(x)
x <- conumee2::CNV.segment(x)

#generate frequency plot
conumee2::CNV.summaryplot(x)
```

## **License**
This package is published under the GPL-2.0 open-source license.

## **Reference**
Please cite our manuscript, if you use CCNV for your research:

*In Preparation*
