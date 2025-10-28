# PCB-Aerobic-Bioaugmentation-Study3
This model simulates the concentration of two PCB congeners: PCB 4 and PCB 19 in both the aqueous and gas phases within laboratory-scale bioreactors containing PCB-contaminated sediment. The script incorporates experimental passive sampler measurements and uses a least-squares minimization approach to estimate congener-specific sampling rates for both passive samplers. Additionally, the model can simulate the decrease in volatile PCB release from sediment to air due to the presence of an aerobic PCB-degrading microorganism, with a known biotransformation rate for PCBs. There are also codes to perform ANOVA and Tukey test to individual PCB congenrs, total PCBs and LCPCBs, as well as for cDNA, DNA and transcript gene ratio for the individual sampling time points.

--------------------------
SHARING/ACCESS/ATTRIBUTION LICENSE INFORMATION
--------------------------

Licenses/restrictions: licensed under the 2-Clause BSD License - see the [LICENSE](LICENSE) file for details.

----------------------
General Information
----------------------

Deposit Title: Polychlorinated Biphenyl (PCB) Reactive Transport Model

This README file was generated on March 10, 2025 by Andres Martinez.

Contributor information:

Andres Martinez, PhD
University of Iowa - Department of Civil & Environmental Engineering
Iowa Superfund Research Program (ISRP)
andres-martinez@uiowa.edu
ORCID: 0000-0002-0572-1494

Principal Investigator: Timothy Mattes, PhD
Principal Investigator email: tim-mattes@uiowa.edu

This work was supported by the National Institutes of Environmental Health Sciences (NIEHS) grant #P42ES013661.  The funding sponsor did not have any role in study design; in collection, analysis, and/or interpretation of data; in creation of the dataset; and/or in the decision to submit this data for publication or deposit it in a repository.

This R project is part of the paper: XXX

Subject: Polychlorinated Biphenyls; Contaminant fate and transport; Paraburkholderia xenovorans LB400; Kinetic phase passive sampling; Bioremediation; Biodegradation; Biosurfactants; Bioavailability; GC-MS/MS

GeoLocation: The sediment used in this study was taken from New Bedford Harbor, MA. Laboratory and analytical work was done at the University of Iowa in Iowa City, IA, USA.

--------
PREREQUISITES & DEPENDENCIES
--------
This section of the ReadMe file lists the necessary software required to run codes in "R".

Software:
- Any web browser (e.g., Google Chrome, Microsoft Edge, Mozilla Firefox, etc.)
- R-studio for easily viewing, editing, and executing "R" code as a regular "R script" file:
https://www.rstudio.com/products/rstudio/download/

--------
SOFTWARE INSTALLATION
--------

This section of the ReadMe file provides short instructions on how to download and install "R Studio".  "R Studio" is an open source (no product license required) integrated development environment (IDE) for "R" and completely free to use.  To install "R Studio" follow the instructions below:

1. Visit the following web address: https://www.rstudio.com/products/rstudio/download/
2. Click the "download" button beneath RStudio Desktop
3. Click the button beneath "Download RStudio Desktop".  This will download the correct installation file based on the operating system detected.
4. Run the installation file and follow on-screen instructions. 

--------
R FILES AND STRUCTURE
--------
It is recommended to create a project in R (e.g., PCB-Aerobic-Bioaugmentation-Study3.Rproj). Download the project file (.Rproj) and the R subfolder where the scripts are located, and the Subfolders.R file. Run first the Subfolder.R file, which will generate all the subfolders for this project.
The structure of this project includes an R subfolder where all the R scripts are located, as previously indicated. There is a Data subfolder where the data are storage, and then an Output subfolder, where the results are located.

--------
Data
--------

SPME and PUF data for individual PCB4 and PCB19 can be found at https://doi.org/XXX. The file 06_Dataset_final_PCBmass.csv needs to be downloaded and saved in the Data folder of the project.


