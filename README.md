# Readme

This repository holds the code for the paper "A practitioners guide to estimating the instantaneous reproduction number, Rt" by Chad W. Milando, Kaitlyn Johnson, Christine Sangphet, Sutyajeet Soneja, Pragati Prasad, Md. Sakhawat Hossain, Benjamin Singer, Harry Hochheiser, and Laura F. White

## List of files

To find packages, we used these scripts:

- `GitHub_search.py`: Python code to search github repositories
- `R_packagesearch.R`: R code to search CRAN and r-universe.dev

Then here is the code for the images in the manuscript:

Figure 3:

- `01_ReportsInfections.R`
- `02_FixedSlidingWindows.R`
- `03_RandomWalk.R`
- `04_Filtering.R`
- `05_B-splines.R`
- `06_GaussianProcess.R`

To run all

- `07_mergePlots.R`

Then to make Supplemental Figure S1:

- `08_ern.R`
- `09_EpiInvert.R`
- `10_estimateR.R`
- `11_mergeSFig1.R`

Then to make the supplemental GP figure:

- `12_GP_demo.R`
