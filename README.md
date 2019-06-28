# diversity-productivity-reforestation
Code used to generate results, figures and tables for DOI: 10.1111/geb.12962

This project uses **forest inventory data** and **plant functional trait data**, processing them in **R scripts** to generate results and figures.

**Forest inventory data**
The forest inventory data used to generate this publication were obtained from CSIRO, Greening Australia, Greenfleet and the Queensland Queensland Department of Agriculture and Fisheries. The CSIRO data is publically accessible (available through the Australian Government Data Access Portal and National Biomass Library). The rest of the data were obtained via data agreements and cannot be reproduced here.

I have provided anonymous plot-level aggregate data, the minimum required to reproduce all results and figures used in this project. I have not provided any plant-level data, site coordinates or identifying information.

**Functional trait data**
Functional trait data was obtained from the TRY Plant Trait Database, field collection made by the authors of this paper, field collections made from colleagues, and Australian herbarium records (for maximum height only). Citations to trait databases and publications for TRY functional trait data are made in the Methods section of the paper. Raw traits from TRY are not reproduced here.

Field collected trait values for specific leaf area and wood density made by Timothy Staples are available here, and will be available via other plant trait databases in the future. Details of field sampling protocol are available via the Method sections of the publication associated with this repository.

Trait values are available as .csv files containing species means. They are located in this repository under the "Functional trait data" subdirectory. These traits are made freely available for future research, but we request that the publication attached to this repository be cited as the source of these values if some or all of them are used:

*Staples TL, Dwyer JM, England JR and Mayfield MM (2019). Productivity does not correlate with species and functional diversity in Australian reforestation plantings across a wide climate gradient. Global Ecology and Biogeography. doi: 10.1111/geb.12962.*

**R scripts** Two scripts are provided in this repository, along with summarised data used for for the statistical models and analyses in the manuscript. These scripts makes extensive use of the code folding functionality in R-studio (Alt + O is the default shortcut on Window and Linux machines to collapse all folds).

**reproduce-results.R** is a minimum working script that contains all of the data and analyses to produce the figures and tables found in the manuscript. It does not, in most cases, show how the raw data was transformed from the raw functional trait and forest inventory data. It makes use of files in the "Data" and "Functions" sub-folder, so please point your working directory to the location of this script and include all other directories in the repository.

**full-processing-and-analysis.R** shows the complete data analysis pathway, but will require the
raw forest inventory and functional trait data acquired and aggregated by the authors of this publication. These are not provided here. This script is provided for transparency and reproducibility, that data processing steps are fully documented.
