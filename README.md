# Beyond adherence thresholds: A simulation study of the optimal classification of longitudinal adherence trajectories from medication refill histories

## Summary

This repository contains the `R` code and report on a simulation study of a sliding-window approach to analyze longitudinal adherence trajectories from electronic healthcare data (EHD). 
The code is released under a [GPL v3](https://www.gnu.org/licenses/gpl-3.0.en.html) license.

The aims of this project were to analyze the effect of different adherence estimation methods, sliding window parameters, and sample characteristics on the performance of a longitudinal clustering algorithm to

a) identify temporal adherence patterns from EHD and

b) classify individuals accurately into the identified clusters. 

The objective of this study is to provide guidance for future longitudinal adherence studies using EHD.

## Contents of this repository

The **accompanying report** (in the `./Report/` directory) describes in detail the simulation study and its results.
The report is the pre-print version of the published version.

The **`R` script** (in the `./Code/` directory) contains the code used for the simulation study, including the custom functions to simulate refill histories.
If you would like to perform a similar simulation study, you may use and adapt the provided script.

If you use (parts of) the `R` scripts, please cite this work and include links to this repository: (https://github.com/Masswear/BeyondThresholds) 

## Feedback

We love your feedback! Although we did our best to comment the `R` code and test it on various systems, it has not been created with an application beyond this work in mind. 
Any comments, suggestions or bug reports are welcome, either through GitHub's own issue reporting options or by e-mail to <s.allemann@unibas.ch>. 

## Thank you

Samuel Allemann
Lyon, France, August, 2018


