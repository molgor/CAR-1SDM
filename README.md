# Conditional autoregressive models for single species distributions using presence-only data
This repository implements the presence-only models developed in the paper: 
[[doi:https://doi.org/10.1101/2021.06.28.450233 ][A joint distribution framework to improve presence-only species distribution models by exploiting opportunistic surveys]]

The folder 'notebooks' contains the interactive pipelines for fitting and displaying models I, II and III. In the current version, the model fitting is performed in R with the package[[https://cran.r-project.org/web/packages/CARBayes/index.html][CARBayes]] (Lee, 2015) and the methods for displaying the results are done in Python. Therefore, it is recommended to use Jupyter notebooks with the R kernel installed for maximum compatibility.

This repository has two purposes. One, to create a reproducible environment for generating the results exposed in the paper (section: application). Second, to serve as a /user's guide/ for adapting the models to different types of problems. 

As this repository is only for demonstrating purposes, it is not intended to be continuously developed in the long term. Must of the ideas specified in models I, II and III are already implemented in the package [[https://github.com/molgor/JSDM4POD][JSDM4POD]] (i.e. /joint species distribution models for presence-only data/). Currently: models I and II are contained in the JSDM4POD. In the coming future model III will be added.

Thank you for your interest in this project. 

Please consider sending your feedback to me: j.escamillamolgora[-at-]lancaster.ac.uk
