# Joint Modeling with a Defective Gompertz Distribution and Frailty

## Overview

This folder contains the R scripts, data, simulation routines, and results associated with the development and evaluation of joint models for longitudinal and survival data using a defective Gompertz distribution with frailty.

The project is organized into separate directories to facilitate code maintenance, reproducibility, and the distinction between simulation studies and real-data applications.

## Directory Structure

### `application/`

Contains the scripts used for the real-data application, including data processing, model fitting, statistical inference, and the generation of tables and figures. These scripts are intended to demonstrate the application of the proposed modeling framework to empirical data.

### `simulation/`

Contains the scripts used to conduct Monte Carlo simulation studies. These files implement the simulation scenarios, generate synthetic datasets, fit the proposed models, and evaluate the resulting parameter estimates and model performance measures.

### `data/`

Contains the datasets required for the real-data application and, when applicable, input files used in the simulation studies. Data provenance and any applicable access or confidentiality restrictions should be documented where relevant.

### `results/`

Contains selected outputs from the application and simulation studies, such as summary tables, figures, parameter estimates, and model evaluation results. Large intermediate files and temporary outputs may be excluded to keep the repository manageable.

## Software and Requirements

The analyses are primarily implemented in R. Depending on the scripts being executed, additional R packages and model-fitting tools may be required. The specific dependencies and their versions should be documented to facilitate reproducibility.

## Reproducibility

The scripts should be run from the appropriate project directory, ensuring that the required input datasets and dependencies are available. Simulation settings, random seeds, and model specifications should be recorded whenever applicable.

Please ensure that all datasets can legally be shared before uploading them to a public repository.

