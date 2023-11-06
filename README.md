# qDATA
### R based analysis of qRT-PCR data

<p align="justify">

### **qDATA** is a bioinformatics tool created with R Shiny in the Drosophila Laboratory, Department of Genetics, Faculty of Biology, University of Bucharest.

### For citing our tool, please use: TBP

#
#
#### Dependencies:

Prior to data analysis we highly recommend installing the latest versions of **R** (https://cran.r-project.org/) and **RStudio** (https://posit.co/downloads/). Our application was designed to accommodate a diversity of operating systems assuming that the R version 4.1.0 is compatible. For installation instructions follow the attached [R_installation.pdf](https://github.com/A-Ionascu/qDATA/blob/main/R_installation.pdf) file.

#

#### Running **qDATA**:

**qDATA** does not require any additional dependencies installation rather than R and RStudio. For Ubuntu-based systems, run [Linux_setup.sh](https://github.com/A-Ionascu/qDATA/blob/main/Linux_setup.sh) prior to running the software. We recommend opening the [qDATA.R](https://github.com/A-Ionascu/qDATA/blob/main/qDATA.R) file in RStudio and running the application with the `Run App` button. Internet connection is required only when the application loads for the first time as all required R packages are downloaded on the local machine.

#

#### Using **qDATA**:

When the application loads, the user is greeted with a modern GUI structured on two windows. The former is the main data analysis environment and the latter offers information about all the abbreviations used in the app. The main window contains an always-on left side panel and a central area for results output. The left side panel contains a browse button which is used to select the input table from any directory on the local machine. In order to change the input data table, the user simply uploads another table using the `Browse` button and the application retrieves the output for the newly uploaded data. Throughout the data analysis process, the user is able to modify parameters from the left side panel with impact on the entire analysis. When an `Update results` button appears under the parameter, the modification is applied after clicking the button. 


**qDATA** can accommodate any number of genes of interest paired with any number of biological or technical replicates. Curently, our implementation is designed for gene expression experiments that rely on a single housekeeping (HK) gene. The user must create a **.csv file** containing the raw data similar to the [input_model.csv](https://github.com/A-Ionascu/qDATA/blob/main/input_model.csv) example. For the first four columns, consistency is required as the script is case and blank space sensitive. The `BR` and `TR` columns are associated with biological replicate and technical replicate numbers. The `Ct` column can be imported from the PCR output table with no limit for the number of decimals.


The qDATA software can be run in the browser by clicking the `Open in Browser` button at the top of the RStudio IDE window. Therefore, the user can use multiple independent instances of the application in its browser mode. When the application is run from the RStudio IDE, the local machine works as the server, both for the RStudio window and the browser instances. Closing the RStudio environment automatically disconnects the server and freezes the browser instances. In order to close the RStudio qDATA window, we recommend using the `Stop` button located at the top of the RStudio IDE Console rather than simply closing the running GUI.


#
#### Contact and feedback:

For offering feedback or any type of inquires about the application, please contact us at **aionascu.g@gmail.com**.  


</p>
