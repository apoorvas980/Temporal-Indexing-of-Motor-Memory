**This folder contains the analysis-ready datasets used by the MATLAB scripts in the main repository.**

1. **behavioral_explicit_lme.csv**

Analysis-ready training data for the externally-cued experiment.

Used by context_analysis_lme.m.

Main variables:

Subject: anonymized participant identifier

HA: hand angle

Time: trial progression variable used in the regression model

Rot_t: current temporal context predictor

Rot_t1: previous-trial perturbation predictor

Rot_t2: two-trials-back perturbation predictor

2. **behavioral_implicit_lme.csv**

Analysis-ready training data for the internally-monitored experiment.

Used by context_analysis_lme.m.

The variables are the same as those in behavioral_explicit_lme.csv.

3. **internallymonitored_model.csv**

Analysis-ready trial-level data for the internally-monitored experiment.

Used by timing_models.m.

Variables:

id: anonymized participant identifier

subTrial: trial number within participant

phase: experimental phase

goDelay: movement preparation interval in milliseconds

rotDir: perturbation direction

handAng: hand angle in degrees

rtGoSig: reaction time relative to the go signal

4. **externallycued_model.csv**

Analysis-ready trial-level data for the externally-cued experiment.

Used by timing_models.m.

Variables are same as in internallymonitored_model.csv

The modelling script performs the required baseline correction and sign alignment internally before fitting the models.
