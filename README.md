# Observation errors in fine-scale DND data
Data and code on observation errors in fine-scale detection/non-detection (DND) data in plant population monitoring.
This work was led by Amaya Richer, Jan Perret, Aurélien Besnard, Anne Charpentier and Guillaume Papuga, in the perspective of being published in a scientific journal. 

# Author contributions
Amaya Richer, Jan Perret, Aurélien Besnard, Anne Charpentier and Guillaume Papuga conceived the study and designed the methodology; Jan Perret and Guillaume Papuga organised and conducted the field sessions from which data are issued; Amaya Richer conducted the data analysis and led the writing of the manuscript. All authors critically contributed to the drafts and gave their final approval for publication. 

Data used in this work were collected for and published in the following paper:
Perret, J., Besnard, A., Charpentier, A., & Papuga, G. (2023). Plants stand still but hide: Imperfect and heterogeneous detection is the rule when counting plants. Journal of Ecology. https://doi.org/10.1111/1365-2745.14110

# Running the code
Download all files and folders in a new folder, then run the *DND_data_Analysis_and_Figures.qmd* script. This script will load the data from the *data* folder, save formatted data in the *saved_tables* folder, fit the models and print the main figures.

Then, run the *DND_data_Simulation.R* script, which contains the code to simulate fine-scale DND data, run occupancy models and print main figures. Because running the code for simulations and model fitting takes a while, outputs of simulations are already saved and stored in simulation_outputs, if one does not desire to rerun simulation analysis.

Supplementary material is in the *DND_data_SuppMat.qmd* script. Knit the script to a PDF format to visualise supplementary figures.

The code was written using R version 4.x.x, and the versions of the packages I used are available in the file package_versions.txt.

# Contact
amaya.richer@cirad.fr or amayaricher@gmail.com
