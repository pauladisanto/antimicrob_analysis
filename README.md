# Antimicrobial analysis

Set of scripts for analysing antimicrobial activity of multiple compounds. The scripts include analysis of: 

#### Synergism:

This script collects the data from the csv file obtaind from the plate reader, prepares the imput data to use the  ``synergyfinder`` library and makes the interaction analysis between two compounds. 
  An example plot of percentage of inhibition

<img src="demo_images/Synergism_inhibition.png" width="500"/>

  An example plot of interaction between two compounds

<img src="demo_images/Synergism.png" width="400"/>

--------------------------------------------------------------

#### Cytotoxicity

There are two scripts for analysing cytotoxicity. Both can collect the data from the csv file obtaind from the plate reader and prepares the imput data. 
One script will let you to calculate EC50 and EC90 values and also make plots for titration curves using the ``drc`` library.
<img src="demo_images/.png" width="400"/>

The other will let you plot all the data as mean and SD. 

<img src="demo_images/Titration_BF.png" width="400"/>

--------------------------------------------------------------

#### Titration curves in bacteria (planctonic and biofilm)
The analysis is similar to the cytotoxicity, but since the data from the csv file obtaind from the plate reader is different the script prepares the imput data to generate the imput for the ``drc``  library.

  An example plot of a titration curve

<img src="demo_images/Titration_BF.png" width="400"/>

--------------------------------------------------------------

#### Growth curves under the effect of antimicrobials, disruption of biofilms, etc.  

  An example plot of biofilm

<img src="demo_images/CV_biofilm.png" width="400"/>

