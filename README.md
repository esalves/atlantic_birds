# Shrinking Body Size in Atlantic Forest Birds as a Response to Climate Change

This repository contains the data and code for the analysis of body size changes in Atlantic Forest passerine birds in response to climate change. The study investigates the relationship between temperature, food availability, and bird morphology over a 28-year period.

## Project Overview

The core of this project is an R Markdown (`.Rmd`) file, `Analysis/atlantic_birds_ms.Rmd`, which performs a comprehensive analysis of a dataset of Atlantic Forest bird measurements. The analysis includes data wrangling, exploratory data analysis, and Bayesian multi-level modeling to test the hypothesis that rising temperatures and declining food resources are leading to a decrease in the body size of passerine birds.

## Repository Structure

- `Analysis/`: This directory contains the main analysis script, data files, and model outputs.
  - `atlantic_birds_ms.Rmd`: The R Markdown script containing the full analysis.
  - `functions.R`: A file containing helper functions used in the analysis.
  - `ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv`: The primary dataset of bird traits.
  - `BirdFuncDat.txt`: A dataset of bird functional traits (Elton traits).
  - `*.rda`, `*.rds`: Saved R data files, including model objects and processed data.
- `RegisteredReport/`: Contains a one-page summary of the registered report for this study.
- `README.md`: This file.

## Setup and Usage

To reproduce the analysis, you will need to have R and RStudio installed on your system.

### 1. Install R and RStudio

- [Download and install R](https://www.r-project.org/)
- [Download and install RStudio Desktop](https://www.rstudio.com/products/rstudio/download/)

### 2. Install Required R Packages

Open the `Analysis/atlantic_birds_ms.Rmd` file in RStudio. The script will prompt you to install the necessary packages if they are not already installed. The required packages are:

- `tidyverse`
- `ggthemes`
- `brms`
- `lme4`
- `ape`
- `geiger`
- `MCMCglmm`
- `mice`
- `tidybayes`
- `reshape2`
- `data.table`
- `metafor`
- `ggmap`

You can install these packages by running the following command in the R console:

```R
install.packages(c("tidyverse", "ggthemes", "brms", "lme4", "ape", "geiger", "MCMCglmm", "mice", "tidybayes", "reshape2", "data.table", "metafor", "ggmap"))
```

### 3. Run the Analysis

1.  Open `Analysis/atlantic_birds_ms.Rmd` in RStudio.
2.  To run the entire analysis and generate the HTML report, click the "Knit" button in the RStudio toolbar.

**Note:** Some of the Bayesian models in the analysis are computationally intensive and may take a long time to run. The pre-computed model objects are saved as `.rda` files in the `Analysis/` directory to allow for faster exploration of the results. By default, the Rmd file is set to use these pre-computed models (`eval=FALSE` on the model fitting chunks).

## Data

- **`ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv`**: This dataset contains morphological measurements, collection dates, and location data for Atlantic Forest birds.
- **`BirdFuncDat.txt`**: This file contains functional trait data for birds, including diet information, from the EltonTraits database.

## Citation

If you use this code or data in your research, please cite the original study. The citation information can be found in the generated report.