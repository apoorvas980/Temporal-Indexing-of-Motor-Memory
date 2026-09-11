**Temporal Indexing of Motor Memory**

This repository contains analysis code and analysis-ready datasets for the manuscript “Temporal Indexing of Motor Memory” by Apoorva Sharma, Hanna Hillman, and Samuel D. McDougle.

The two MATLAB analysis scripts are located in the main repository folder, and all datasets are located in the data folder.

**Files in the main folder:**

1. **context_analysis_lme.m**
Runs the mixed-effects models and participant-level regression analyses for the internally-monitored and externally-cued experiments.

2. **timing_models.m**
Runs the computational modelling analyses for the internally-monitored and externally-cued experiments. The script fits the Subjective Timing Model and the Granule-Cell Basis Model, evaluates test-phase generalization, and reports the model comparison statistics used in the manuscript.

_Data files used by context_analysis_lme.m:_
a. data/behavioral_explicit_lme.csv
b. data/behavioral_implicit_lme.csv

_Data files used by timing_models.m:_
a. data/internallymonitored_model.csv
b. data/externallycued_model.csv

How to run:

1. Open MATLAB in the main repository folder.

2. Run:
context_analysis_lme to reproduce the behavioral mixed-effects and participant-level regression analyses.

3. Run:
timing_models to reproduce the computational modelling analyses.

Both scripts automatically load the required files from the data folder.

**Software requirements:**
MATLAB R2023a
Statistics and Machine Learning Toolbox
Optimization Toolbox

The scripts print the relevant statistical results to the MATLAB Command Window and do not save additional output files.

See the README inside the data folder for descriptions of the datasets and variables.
