## Macroecological factors shaping the geographic distribution of *Spirodela polyrhiza*

This repository serves to store R scripts to reproduce analyses and figures presented in the manuscript "Macroecological factors shaping the geographic distribution of *Spirodela polyrhiza*, with reflections on predictor selection in ecological niche modeling".

Authors: Marlon E. Cobos*, A. Townsend Peterson

Contact: manubio13@gmail.com


### Description
This repository contains scripts to:

- Download relevant data.
    - Occurrence data from GBIF and BIEN.
    - Environmental data from WorldClim, SoilGrids, and Global Gridded Soil Phosphorus.
- Prepare datasets to run analyses.
    - Clean and filter occurrence data.
    - Generate areas for model calibration.
    - Process and generate raster environmental variables.
    - Prepare data for initial analyses and ecological niche modeling.
- Run analyses for the project.
    - Graphical explorations of data.
    - Correlation analyses.
    - Model calibration (model selection).
    - Final models and transfers.
- Create figures to interpret results.
    - Figures presented in the manuscript.
    - Other figures found in the supplementary materials.


### Scripts

- <a href="https://github.com/marlonecobos/Spirodelap_ENM/blob/main/Scripts/00_Functions.R" target="_blank">Script</a> containing functions.
- <a href="https://github.com/marlonecobos/Spirodelap_ENM/blob/main/Scripts/Complete_process.R" target="_blank">Script</a> to run the complete set of analyses needed to reproduce results from the manuscript.
