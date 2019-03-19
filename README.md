# Origin-of-a-preferential-avulsion-node-on-lowland-river-deltas
Code of the numerical model presented in "Origin of a preferential avulsion node on lowland river deltas" by Chadwick et al., (2019)

Austin Chadwick
3/19/19
achadwick@caltech.edu

This folder contains code for the manuscript "Origin of a preferential avulsion node on lowland river deltas." Data was produced by running this code repeatedly, changing the input parameters as described in the manuscript. Below is a brief outline of what each file does in this folder.

MultiDimless.m 
This is the primary script, responsible for running repeated avulsion cycles under a given set of input conditions. (Title is shorthand for "Multiple avulsions, dimensionless morphodynamics")

Dimless.m
This is the secondary script, responsible for morphodynamics for an individual avulsion cycle. This is called by "MultiDimless" script (Title is shorthand for "dimensionless morphodynamics")

interp1qr.m, intersections.m, newton.m, polyfitZero.m, randp.m
These matlab functions were download from file exchange online, and are used by Dimless.m and MultiDimless.m. Austin did not write this code. 
