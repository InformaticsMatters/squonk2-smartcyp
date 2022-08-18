This describes how to run the `  smartcyp-prediction` job from the `comp chem` category in the `smartcyp` collection.

## What the job does

This job calculates predictions for sites of cytochrome p450 metabolism in molecules.

SMARTCyp is a tool from the University of Copenhagen that predicts the atoms that are metabolised by Cytochrome P450
enzymes. It can be used as a guide to whether your molecules might be rapidly metabolised.

The SMARTCyp calculation is slower than many other molecular predictions, and most suited to small numbers
(100’s or a few 1000's) of molecules.

## Implementation details

* Job implementation: [Predictor.java](/app/src/main/java/squonk/jobs/smartcyp/Predictor.java)
* Job definition: `jobs.  smartcyp-prediction` in [rdkit.yaml](../smartcyp.yaml)
* About SmartCYP: https://smartcyp.sund.ku.dk/

Our implementation is based on the Java source code provided by the University of Copenhagen. We have adapted it slightly
to function as a Squonk job, but there are no material changes to the algorithm and the results should be identical
to those generated by the SMARTCyp web site or download versions.

One enhancement we do make is related to structures with multiple fragments. The University of Copenhagen implementation
breaks the structure into fragments and keeps the largest fragment and outputs results for just that largest fragment.
The Squonk implementation also calculates scores only for the largest fragment, but maps those scores back to the original
structure, thus not changing the original structure in any way.

SMARTCyp uses an old version of CDK (version 1.4.8) in its implementation. Squonk's SMARTCyp implementation also uses
this version to ensure reproducibility of results. How this means that there is a small possibility that molecules
handled by other CDK based jobs (using a newer version of CDK) may not be handled correctly by this job
(e.g. if there are bugs or changes in how molecules are represented between the CDK versions).

## How to run the job

### Inputs

* **Input molecules**: molecules to predict.

### Options

* **Output file** - name for the output SD-file
* **Calculate general reactivity** - calculate the general reactivity
* **Calculate 2D6 reactivity** - calculate the 2D6 reactivity
* **Calculate 2C9 reactivity** - calculate the 2C9 reactivity
* **Apply empirical N-Oxidation corrections** - apply empirical N-Oxidation corrections to calculations
* **Maximum rank** - output only this number of predicted sites for each molecule
* ** Score threshold** - don't output sites below this score

## Output

The following properties are added to the output SD-file, if those predictions are selected:

* SMARTCyp_GEN - general P450 metabolism (e.g. 3A4)
* SMARTCyp_2D6 - prediction for the 2D6 isoform
* SMARTCyp_2C9 - prediction for the 2C9 isoform