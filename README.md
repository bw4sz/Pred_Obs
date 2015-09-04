
#Readme file for 'Reduced co-occurrence among closely related hummingbird species'	

MS submitted as a Research Paper to Global Ecology and Biogeography 

Ben G. Weinstein, corresponding author, email: benjamin.weinstein@stonybrook.edu
Dept. of Ecology and Evolution, Stony Brook University, Stony Brook, New York 11794; USA

Juan Luis Parra, email: juanl.parra@udea.edu.co 		
Grupo de Ecología y Evolucíon de Vertebrados, Instituto de Biología, Facultad de Ciencias Exactas y Naturales, Universidad de Antioquia, Medellín, Colombia

Catherine H. Graham, email: Catherine.Graham@stonybrook.edu
Dept. of Ecology and Evolution, Stony Brook University, Stony Brook, New York 11794; USA

For the greatest transparency and reproducibility, we have included all R files used to create this analysis. The entire repository can be found here: https://github.com/bw4sz/Pred_Obs. For git users, the repo can be cloned

```{}
git clone https://github.com/bw4sz/Pred_Obs.git
```

The supplamental documents are labeled by their filename and their order. On the Journal website, they are just listed as Appendix 1,2,3 etc, but follow the order below.

----------

* Appendix 1: SpeciesOverlap.html

This is the main section of analysis which takes in the species assemblages, localities and runs the niche models (SDMUpdated.R). The predicted habitat suitability maps are used to create predicted environment assemblages at each geographic location for which we have observed assemblage lists. The phylogenetic distance to the closest related species is calculated for each species in each assemblage and parameters are fit using hierarchical bayesian approaches. The posterior estimates are compared to simulated posteriors created in Appendix 2

* Apppendix 2: SimulationsPezEcho.html

Using the scape function from the R package PEZ, we simulate different occurrence patterns based on models of trait evolution that included phylogenetic signal and/or repulsion in an occurrence trait. We test our result for sensitivity to tree size and shape.

* Appendix 3: SDMEval.html

Evaluation metrics of our ensemble niche models for all species, as well as correlation values among test statistics and modeling methods. We also compared the sensitivity and specificity of species prediction based on thresholds of habitat suitability and found a 0.05 quantile best minimized overprediction and underprediction. 

* Appendix 4: Sympatry.html

Analysis of the severity of spatial patterns of geographic barriers between sister taxa and time since divergence. We also fit mixed effects models to explore the statistical relevance of this relationship at multiple phylogenetic depths.

* Appendix 5: Non_Independence.html

Fitting alternative bayesian, glm, and glmmm models to the data and comparing results to randomization assemblages based on tip swapping. These approaches evaluate whether we could arrive at our result based solely on the topology and relatedness of the taxa.

* Appendix 6: README.html

ReadME file you are currently reading orienting the reader to the goals and relavance of each of the appendices. 

