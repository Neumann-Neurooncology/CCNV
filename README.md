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
devtools::install_github("AntoniaGocke/CCNV")
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

## **License**
This package is published under the GPL-2.0 open-source license.

## **Reference**
Please cite our manuscript, if you use CCNV for your research:

*In Preparation*
