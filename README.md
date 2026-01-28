---
title: "Documentation of DILI Scoring System"
author: Haochen Li
date: "2026-01-07"
---


# Overview

## Background and Motivation 
Prior to this work, DILI assessment relied heavily on subjective interpretation, particularly through manual inspection of data visualizations such as histograms. This approach has several limitations:

1. High subjectivity\
Histogram construction (e.g., binning strategy) is inherently subjective and sensitive to experimental variance, which can substantially influence interpretation.

2. Limited comparability across compounds\
When toxicity differences between compounds are subtle, it becomes difficult to make clear, quantitative comparisons using individual assay readouts or visual inspection alone.

3. Over-reliance on individual concentration levels\
In conventional analyses, concentration levels that do not reach statistical significance are often excluded, unless a strong overall dose-dependent trend is observed. As a result, potentially informative signals may be overlooked.

## Objective

This code is designed to construct and evaluate a DILI scoring system based on organchip experimental data, with the goal of providing a quantitative and comparable assessment of compound-level hepatotoxic risk.

Specifically, the system takes multiple biological assay readouts as inputs—including but not limited to cell viability, albumin, ATP, and enzymatic markers—measured across different treatment conditions and concentration levels. These data are systematically integrated, normalized, and transformed to produce interpretable compound-level DILI scores.

The resulting scores can be used to:

- Rank and compare hepatotoxic risk across compounds;

- Assess toxicity differences under different treatment conditions or experimental setups;

- Evaluate consistency and reproducibility across experimental batches or laboratories.


## Run the code
In MAC -> Terminal-> Run  "R --vanilla"

After Entering R: 

setwd("/yourpath/Scoring System Program")
p'l
source("run.R")

Follow the interactive prompts to generate toxicity scores for each compound and assay at the selected time points. 
The results will be saved as Excel files.

After the scoring results are generated and organized in a wide format, run:

source("heatmap.R")

Select the input file and columns interactively, adjust color limits and scaling if needed, 
and the program will generate a toxicity heatmap image saved to the output/ directory.
