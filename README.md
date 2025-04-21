# StateSpaceMisspecification (ssm.R)


**phylogenetic simulation and reconstruction comparison pipeline**

This repository provides a pipeline for simulating phylogenetic trees and associated datasets, reconstructing trees under various models (partitioned/unpartitioned, likelihood/parsimony), and comparing their performance.

---

## 🔧 Features

- The pipeline:
  - Simulate pure birth or birth-death phylogenies with [pbtree](https://www.rdocumentation.org/packages/phytools/versions/2.4-4/topics/pbtree)
  - Generate Mk-model character matrices with [sim.Mk](https://www.rdocumentation.org/packages/phytools/versions/1.9-16/topics/sim.history)
  - Reconstruct phylogenies using [RAxML](https://github.com/stamatak/standard-RAxML) (and optionally [PAUP*](https://paup.phylosolutions.com/))
  - Compare reconstructed vs. true trees with RF, PID, and CID distances with [TreeDist](https://ms609.github.io/TreeDist/)
- Parallel execution support

---

## 📦 Requirements

Install the following R packages before use:

```r
install.packages(c(
  "ape", "phytools", "TreeDist", "phangorn",
  "processx", "foreach", "doParallel", "Matrix", "readr"
))
```

[RAxML](https://github.com/stamatak/standard-RAxML) (and [PAUP*](https://paup.phylosolutions.com/) if intend to include parsimony) must also be installed and accessible via command line.

---

## 🚀 Example

See [`run_example.Rmd`](examples/run_example.Rmd) for an example call with parameters explained.

See [`data_sim.Rmd`](examples/data_sim.Rmd) for templates to perform and summarize simulation

---

## 📊 Visualization

Plots for comparing results (e.g., `plot_ssmr`, `plot_ssmd`) are provided in [`run_example.Rmd`](examples/plot_ssm.Rmd).

---


## 📄 Citation

If you use this package, please cite:

> (coming soon...)

---
