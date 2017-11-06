# qPCR-VA
Program and supplements to visually analyze qPCR-data
Content:
1. Requirements
2. Install
3. Usage

###########################################################################################

1. Requirements before Installation:

- OS: Windows 64-Bit

- R version 3.4.2 or higher(get newest version at:https://www.r-project.org/)

- optionaly download R-Studio at: https://www.rstudio.com/products/rstudio/download/	(only if you intend to change the R-script working in the background)

- jre version 1.8 or higher(get newest version at:http://www.oracle.com/technetwork/java/javase/downloads/index.html)


#############################################################################################

2. Install:

- Install jri and rJava by opening R or R-Studio and typing "install.packages('rJava')"

- Install rmngb with "install.packages('rmngb')"

- Install reshape2 with "install.packages('reshape2')"

- Install ddCt with "source('https://bioconductor.org/biocLite.R')"
		    "biocLite('ddCt')"
			
			
####################################      ERROR HELP      ###################################
If you get an Error saying something like: "could not move temporary file from 'C:\\Path' to 'C:\\Path'"
try typing: debug(utils:::unpackPkgZip)	
and then installing the package again --> You will have to press Enter a lot until it should say: "package downloded
succesfully ...something... check sum MD5"
Sometimes it can be necessary to repeat the step mentioned above to get rid of any issues					

Check if all necessary packages are downloaded correctly, by typing:	
	library("rmngb")
	library("reshape2")	
	library("ddCt")	
	library("ggplot2")	
If there is a Warning like could not find package "NameOfPackage" try installing it manually with "install.packages('NameOfPackage')"
#############################################################################################

- Execute Install.bat by double clicking and type in the according paths. The necessary libraries will be copied to your JRE home directory

- Run PlotAndStats.jar

#############################################################################################

3. Usage:

- This program is intended to work with data generated by LinReg. The Resulting Excelsheet should contain 4 Header rows with general Information and then the Data.


- Please make sure the Colums are in the following order:
  
  name    indiv_PCR_eff    Amplicon    threshold    mean_PCR_eff    cq    N0


- The names must have the format: Wellplate_Gene_Condition_Group		(e.g. A1_TBP_0Gy_15min)


- After controling your Excelsheet start the Programm and fill in the the mask according to the data of your Excelsheet
  
  *if you used spacebars or other separators (apart from "_") in the Excelsheet also use them while filling in the Mask

- When the whole mask is filled out press 'Apply'. The program wll now create a bar plot and show the coressponding data in a Table underneath.

- If you wish to save the Plot and Table press the save button and select the directory where the files should be saved

#########################################################################################################

Created by:

Justin Murtagh
