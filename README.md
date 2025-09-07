# Lipolook

Lipolook is a primary statistical analysis tool built in R Studio specifically for LC-MS lipidomics data. 
Detailed explainations of pipeline usage can be found in the *Lipolook_instructions.R* document provided in this repository. Open with R Studio and choose the `knit` option on the top ribbon to produce the user friendly HTML file. 
The specified file path should contain two .csv files: *raw data* and *lipids tested*

#### Capabilities

#### Inputs

An instructional HTML file (*Lipolook_instructions.R*) is provided for the data input and cleaning sections of the pipeline. Everything past these sections should be fully automated. However, if any issues do arise, a separate annotated copy of the script (*Lipolook_script_annotated.R*) is given for information on what each section of code does. 

#### Outputs
Outputs will be saved in the following directory structure: 
<img width="2521" height="3123" alt="File structure" src="https://github.com/user-attachments/assets/2c358636-75d7-4d3e-95a2-416dcbc1808b" />









Effort must be made to emphasise that this is a primary version of this pipline. While every effort has been made to automate it as well as possible, more data is needed to refine its capabilities. Therefore, when using the pipeline, please do flag any issues found to sullivan-al-kadhomiyj@cardiff.ac.uk. 

## Version Control

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
