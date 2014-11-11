Reduced co-occurrence among closely related species: Evidence from simulated, predicted and observed hummingbirds assemblages.
----------------------------------------------------

*MS submitted Ecology in the “Article” category.*

Ben G. Weinstein, corresponding author, email: bweinste@life.bio.sunsysb.edu
Department of Ecology and Evolution, Stony Brook University, Stony Brook, New York 11794; USA
Juan Luis Parra, email: juanluisparra@gmail.com
Instituto de Biología, Facultad de Ciencias Exactas y Naturales, Universidad de Antioquia, Medellín, Colombia
Catherine H. Graham, email: Catherine.Graham@stonybrook.edu
Department of Ecology and Evolution, Stony Brook University, Stony Brook, New York 11794; USA

Code
--------

The main function is SpeciesOverlap.R, it is both a wrapper for source functions and a central script to create most of figures and modeling outputs
- Reads in input data
- Calls niche modeling script SDMupdated.R given a cell size
- Retrieves the niche model results and computes distance matrices for all assemblage types
- Computes logistic regression for all assemblage types
- Creates figures and prints to file

Many of the minor functions are in SpeciesOverlapSourceFunctions.R

Simulations are made using simulationfigures.R - a small copy is placed so users may run the code quickly at smale scales. For several hundred iterations, the NSF Stampede was used, and the entire set of files (including the shell scripts for run calls ) is in the folder **cluster**. The source script pglmmrichness.R is mostly a wrapper of the pglmm.sim code from the package *picante.*

As always, please check the relative path names depending on where this git repo is cloned. Many packages are required, see the tops of each script for specifics. 
