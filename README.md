# README #

C-Application, currently running in XCode

### What is this repository for? ###

* This is my master-thesis in computational physics at the chair of Prof. Hartmut Ruhl at LMU-Munich.
* In detail, plasma simulations are implemented, using FDTD algorithm and a Hybrid Field approach to propagate fields on the numerical grid. Particle pushes are performed by using a Boris Pusher or alternatively an adaptive Nyström Pusher, respecting the relative effects on each particle.
* The goal of this thesis is to find a criteria of when to be able to cut off the collected history (x, t) for individual particles. This promises to save storage space and computational power for large scaled simulations.

### Setup ###

* Run "HybridFields.xcodeproj" and check system path in "system(...)" command
* Configuration --> none
* Dependencies --> none (yet)
* Database configuration
