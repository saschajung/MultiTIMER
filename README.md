# MultiTIMER

MultiTIMER is an R package implementing a multi-tissue transcriptional aging
clock trained on ArchS4 RNA-seq data. 

## Where to start

Note that for best perfomance, new data ought to be concatenated with the training data, intersecting on shared genes (since we do not impute missing gene data). This is achieved by reusing the functions present in `R/ArchS4_helpers.R` so that the concatenated matrix has user test data normalized and batch corrected together with training data, including adding to the `batchid` object in the correct end the number of additional batches in your data. We assume most user test data does not have age metadata, but the training should still perform if the concatenated matrix is subsetted to only the original training data for `trainModel` with the full concatenated matrix then passed to `predictAge`.

Because of this difficulty, an alternative was developed that did not have this complication.

The maintained pipeline for **training, validating, and applying** the clock to
new RNA-seq data is now part of the REVIVE platform:

> **https://github.com/saschajung/REVIVE**

REVIVE provides the packaged workflow (frozen SVA → h2o GLM → Horvath age
back-transform), a minimal data package to retrain the clock, and worked
validation examples. New users should start there.

## This repository

These are the original lower-level training helpers the clock is built from:

- `R/ArchS4_helpers.R` ArchS4 sample selection and batch correction. 
- `R/PredictorTraining.R`  predictor (gene) selection and clock training
- `R/Transformation.R`  Presently only includes identity
- `ExampleScript.R`  legacy entry point that shows how to train the clock on original training data only
