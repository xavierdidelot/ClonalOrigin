# ClonalOrigin
Automatically exported from code.google.com/p/clonalorigin

# Introduction #

Bacteria, unlike us, can reproduce on their own. They do however have mechanisms that transfer DNA between organisms, a process more formally known as recombination. The mechanisms by which recombination takes place have been studied extensively in the laboratory but much remains to be understood concerning how, when and where recombination takes place within natural populations of bacteria and how it helps them to adapt to new environments. `ClonalOrigin` performs a comparative analysis of the sequences of a sample of bacterial genomes in order to reconstruct the recombination events that have taken place in their ancestry.

`ClonalOrigin` is described in the following paper:

Didelot X, Lawson D, Darling A, Falush D (2010) Inference of homologous recombination in bacteria using whole genome sequences. Genetics doi:10.1534/genetics.110.120121 http://www.genetics.org/cgi/content/abstract/genetics.110.120121v1

# Usage #

Instructions for how to download and install `ClonalOrigin` are available at:
http://code.google.com/p/clonalorigin/wiki/Download_And_Install_ClonalOrigin

Instructions for how to use `ClonalOrigin` once installed are available at:
http://code.google.com/p/clonalorigin/wiki/FromGenomeAssemblyToRecombination

# Estimating strength of bias in the recombination process #

We define biased recombination in contrast to free recombination where all individuals in the population are equally likely to recombine. There are many factors contributing to recombination being biased rather than free. Laboratory experiments have shown that the recombination process is homology dependent whereby it tends to happen more often between individuals that are less diverged. Furthermore, the geographical and ecological structures observed in many bacterial populations implies a greater opportunity of recombination for pairs of cells that are closely related. Purifying selection may also effectively prevent recombination between distantly related bacteria. All these effects would clearly be hard to disentangle, and here we group them all under the single concept of biased recombination. The strength of this bias is an important factor to take into account in order to understand recombination in bacteria. In particular, this determines how often recombination happens within the diversity of the population under study rather than from other sources.

We are introducing a model for biased recombination which is based on the `ClonalOrigin` model. We use approximate Bayesian computation and whole genome data to infer the rate of bias in the recombination process in bacteria. The user guide and Matlab code can be downloaded from:
http://www.stats.ox.ac.uk/~ansari/BiasedRecV1.tgz
