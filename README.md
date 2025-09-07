# Lipolook

Lipolook is a primary statistical analysis tool built in R Studio specifically for LC-MS lipidomics data.  
The *Lipolook_script.R* file is the one that should be used for execution of the pipeline with your data.  
Test data for this pipeline is provided in the `test_data` directory.  
Detailed explainations of pipeline inputs can be found in the *Lipolook_instructions.R* document provided in this repository. Open with R Studio and choose the `knit` option on the top ribbon to produce the user friendly HTML file. 
  
### Inputs
  
Inputs for the Lipolook pipeline include `raw_data`, `lipids_tested` and `lipid_categories` files. These files must have a specified layout in order for the pipeline to function correctly. An instructional file (*Lipolook_instructions.R*) is provided for the data input and cleaning sections of the pipeline. This will guide and assist with any manual steps included. Everything past these sections should be fully automated. However, if any issues do arise, a separate annotated copy of the script (*Lipolook_script_annotated.R*) is given for information on what each section of code does. 
  
### Outputs
Outputs will be saved in the following directory structure:   
  
<img width="500" height="800" alt="File structure" src="https://github.com/user-attachments/assets/2c358636-75d7-4d3e-95a2-416dcbc1808b" />  
  
This file structure will be built using the pipeline.  
  
**Main numerical outputs for the Lipolook pipeline are:**  
* Group averages across all lipids
* Histograms for each lipid
* Distribution summaries and statistical confimation/denial or normality for each lipid (across all samples and within groups)
* Distribution summaries and statistical confimation/denial or normality for each lipid (across all samples and within groups) using log transformed data
* Mann-Whitney U test significance results
* Kruskal-Wallis H test significance results
* Spearman correlation results
  
**Main visual outputs for the Lipolook pipeline are:**  
* Histograms
* Forest plots
* Bar plots
* Correlation matrix heatmaps
* Volcano-style plots
  
### Future aims
  
While Lipolook is currently in a primative state, its ultimate aim is to provide people outside of the Lipidomics field with a way to easily analyse their data. While it currently does this in a basic form, future improvements should include:  
* Better automation
* A summary HTML file
* A more comprehensive instructional document
* More statistical options
* Futher speciality for lipidomics data
* A more user friendly interface
  
These changes are currently under development and hopfully will be provided in an updated version promptly.  
  
Effort must be made to emphasise that this is a primary version of this pipline. While every effort has been made to automate it as well as possible, more data is needed to refine its capabilities. Therefore, when using the pipeline, please do flag any issues found to sullivan-al-kadhomiyj@cardiff.ac.uk. 
  
### Version Control
  
* `moments` - 0.14.1
* `tidyr` - 1.3.1
* `dplyr` - 1.1.4
* `ggplot2` - 3.5.2
* `stats` - 4.4.1
* `pheatmap` - 1.0.12
* `ggrepel` - 0.9.6
* `readr` - 2.1.5
* `tidyverse` - 2.0.0
* `FSA`
