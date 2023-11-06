# qDATA
R based analysis of qRT-PCR data


#### qDATA is a bioinformatics tool created with R Shiny in the Drosophila Laboratory, Department of Genetics, Faculty of Biology, University of Bucharest.

#### For citing our tool, please use: TBP


Prior to data analysis we highly recommend installing the latest versions of R (https://cran.r-project.org/) and RStudio (https://posit.co/downloads/). Our application was designed to accommodate a diversity of operating systems assuming that the R version 4.1.0 is compatible. For installation instructions follow the attached R_installation.pdf file.





When the application loads, the user is greeted with a modern GUI structured on two windows. The former is the main data analysis environment and the latter offers information about all the abbreviations used in the app. The main window contains an always-on left side panel and a central area for results output. The left side panel contains a browse button which is used to select the input table from any directory on the local machine. In order to change the input data table, the user simply uploads another table using the `Browse` button and the application retrieves the output for the newly uploaded data.



The qDATA software can be run in the browser by clicking the `Open in Browser` button at the top of the RStudio IDE window. Therefore, the user can use multiple independent instances of the application in its browser mode. When the application is run from the RStudio IDE, the local machine works as the server, both for the RStudio window and the browser instances. Closing the RStudio environment automatically disconnects the server and freezes the browser instances. In order to close the RStudio qDATA window, we recommend using the `Stop` button located at the top of the RStudio IDE Console rather than simply closing the running GUI.
