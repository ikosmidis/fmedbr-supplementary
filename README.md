# Supplementary material for “Focused median bias reduction”
[Davide Benussi](https://davidebenussi.github.io/), [Ioannis
Kosmidis](https://www.ikosmidis.com), [Alessandra
Salvan](https://homes.stat.unipd.it/alessandrasalvan/en/home-2/), [Nicola
Sartori](https://homes.stat.unipd.it/nicolasartori/en/home-2/)
October 8, 2026

# Directory structure

The directory `code/` contains the R scripts that reproduce all the
numerical experiments in the manuscript

> Benussi D, Kosmidis I, Salvan A, Sartori N (2026). Focused median bias
> reduction. https://arxiv.org/abs/2606.28597

and the Supplementary Material document
[`fmedbr-supplementary.pdf`](fmedbr-supplementary.pdf).

Running those scripts populates the directory `results/` with the
results of those experiments.

Then, running the script `produce_outputs.R` generates all tables (in
the directory `tables/` in LaTeX format) and graphics (in the directory
`figures/` in PDF format) reported in the main text and the
Supplementary Material document.

# R version and contributed packages

All results are reproducible using R version 4.6.1 (2026-06-24) and the
contributed packages

<table style="width:44%;">
<colgroup>
<col style="width: 26%" />
<col style="width: 18%" />
</colgroup>
<thead>
<tr>
<th>Package</th>
<th>Version</th>
</tr>
</thead>
<tbody>
<tr>
<td>BAMBI</td>
<td>2.3.7</td>
</tr>
<tr>
<td>brbetabinomial</td>
<td>1.0</td>
</tr>
<tr>
<td>brglm2</td>
<td>1.1.0</td>
</tr>
<tr>
<td>detectseparation</td>
<td>0.4.0</td>
</tr>
<tr>
<td>dplyr</td>
<td>1.2.1</td>
</tr>
<tr>
<td>enrichwith</td>
<td>0.6</td>
</tr>
<tr>
<td>focuson</td>
<td>0.4</td>
</tr>
<tr>
<td>future.apply</td>
<td>1.20.2</td>
</tr>
<tr>
<td>ggplot2</td>
<td>4.0.3</td>
</tr>
<tr>
<td>ISLR2</td>
<td>1.3-2</td>
</tr>
<tr>
<td>likelihoodAsy</td>
<td>0.51</td>
</tr>
<tr>
<td>MuMIn</td>
<td>1.48.19</td>
</tr>
<tr>
<td>mvtnorm</td>
<td>1.4-2</td>
</tr>
<tr>
<td>numDeriv</td>
<td>2016.8-1.1</td>
</tr>
<tr>
<td>ordinal</td>
<td>2026.7-26</td>
</tr>
<tr>
<td>patchwork</td>
<td>1.3.2</td>
</tr>
<tr>
<td>pracma</td>
<td>2.4.6</td>
</tr>
<tr>
<td>progressr</td>
<td>1.0.0</td>
</tr>
<tr>
<td>robustbase</td>
<td>0.99-7</td>
</tr>
<tr>
<td>tidyr</td>
<td>1.3.2</td>
</tr>
<tr>
<td>tinytable</td>
<td>0.19.0</td>
</tr>
<tr>
<td>tinytest</td>
<td>1.4.3</td>
</tr>
</tbody>
</table>

The `focuson` R package is available
[here](https://github.com/ikosmidis/focuson).

The `brbetabinomial` R package is available
[here](https://github.com/eulogepagui/brbetabinomial).

# Reproducing the results

## Path

All scripts specify the path to the supplementary material path as
`base_dir`. This is currently set to `.` assuming that the working
directory in R is set to the current `git` repository. If this is not
the case for your setup, you should set `base_dir` appropriately.

## Details

The following table lists the R scripts that need to be executed in
order to reproduce the results. The table also lists the outputs from
each script, and their label if they are shown in the main text or the
Supplementary Material document. Some of the outputs are intermediate
results, so the scripts should be executed in the order shown.

<table style="width:100%;">
<colgroup>
<col style="width: 40%" />
<col style="width: 53%" />
<col style="width: 5%" />
</colgroup>
<thead>
<tr>
<th>Script</th>
<th>Output</th>
<th>Type</th>
</tr>
</thead>
<tbody>
<tr>
<td><a href="code/beta-binomial.R">beta-binomial.R</a></td>
<td><a href="results/beta-binomial.rda">beta-binomial.rda</a></td>
<td></td>
</tr>
<tr>
<td><a href="code/biv-von-Mises-var1.R">biv-von-Mises-var1.R</a></td>
<td><a
href="results/bvmsin-circular-variance.rda">bvmsin-circular-variance.rda</a></td>
<td></td>
</tr>
<tr>
<td><a href="code/fic.R">fic.R</a></td>
<td><a href="results/fic.rda">fic.rda</a></td>
<td></td>
</tr>
<tr>
<td><a
href="code/mahalanobis-distance.R">mahalanobis-distance.R</a></td>
<td><a
href="results/mahalanobis-distance.rda">mahalanobis-distance.rda</a></td>
<td></td>
</tr>
<tr>
<td><a
href="code/mahalanobis-distance-2.R">mahalanobis-distance-2.R</a></td>
<td><a
href="results/mahalanobis-distance-2sample.rda">mahalanobis-distance-2sample.rda</a></td>
<td></td>
</tr>
<tr>
<td><a href="code/marginal-effects.R">marginal-effects.R</a></td>
<td><a href="results/marginal-effects.rda">marginal-effects.rda</a></td>
<td></td>
</tr>
<tr>
<td><a href="code/multiple-mediator.R">multiple-mediator.R</a></td>
<td><a
href="results/multiple-mediator-a.rda">multiple-mediator-a.rda</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a
href="results/multiple-mediator-b.rda">multiple-mediator-b.rda</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a
href="results/multiple-mediator-c.rda">multiple-mediator-c.rda</a></td>
<td></td>
</tr>
<tr>
<td><a href="code/ordinal-superiority.R">ordinal-superiority.R</a></td>
<td><a href="results/osm.rda">osm.rda</a></td>
<td></td>
</tr>
<tr>
<td><a href="code/weibull.R">weibull.R</a></td>
<td><a href="results/weibull.rda">weibull.rda</a></td>
<td></td>
</tr>
<tr>
<td><a href="produce_outputs.R">produce_outputs.R</a></td>
<td><a href="tables/beta-binomial.tex">beta-binomial.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a
href="tables/bvmsin-circular-variance.tex">bvmsin-circular-variance.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a href="tables/fic.tex">fic.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a
href="tables/mahalanobis-distance.tex">mahalanobis-distance.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a
href="tables/mahalanobis-distance-2sample.tex">mahalanobis-distance-2sample.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a
href="tables/multiple-mediator-a1.tex">multiple-mediator-a1.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a
href="tables/multiple-mediator-a4.tex">multiple-mediator-a4.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a
href="tables/multiple-mediator-b1.tex">multiple-mediator-b1.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a
href="tables/multiple-mediator-b4.tex">multiple-mediator-b4.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a
href="tables/multiple-mediator-c1.tex">multiple-mediator-c1.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a
href="tables/multiple-mediator-c1.tex">multiple-mediator-c4.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a href="tables/osm.tex">osm.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a href="tables/weibull1.tex">weibull1.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a href="tables/weibull5.tex">weibull5.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a href="tables/weibull10.tex">weibull10.tex</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td><a href="figures/marginal-effects.pdf">marginal-effects.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td></td>
<td><a
href="figures/marginal-effects-cover.pdf">marginal-effects-cover.pdf</a></td>
<td>Figure</td>
</tr>
</tbody>
</table>
