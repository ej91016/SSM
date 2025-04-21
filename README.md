# SSM (ssm.R)

**_A pipeline for simulating and testing phylogenetic reconstructions with state-space variation_**

---

## 🔧 Features

- The pipeline:
  - Simulate pure birth or birth-death phylogenies with [pbtree](https://www.rdocumentation.org/packages/phytools/versions/2.4-4/topics/pbtree)
  - Generate Mk-model character matrices with [sim.Mk](https://www.rdocumentation.org/packages/phytools/versions/1.9-16/topics/sim.history)
  - Reconstruct phylogenies using [RAxML](https://github.com/stamatak/standard-RAxML) (and optionally [PAUP*](https://paup.phylosolutions.com/)) with different state space assumption 
  - Compare reconstructed vs. true trees with RF, PID, and CID distances with [TreeDist](https://ms609.github.io/TreeDist/)
- Parallel execution support

---

## 📦 Setup

Install the required R packages:

```r
install.packages(c(
  "ape", "phytools", "TreeDist", "phangorn",
  "processx", "foreach", "doParallel", "Matrix", "readr"
))
```

[RAxML](https://github.com/stamatak/standard-RAxML) (and [PAUP*](https://paup.phylosolutions.com/) if intend to include parsimony) must also be installed and accessible via command line.

---

## 🚀 Getting Started

See [`run_example.Rmd`](examples/run_example.Rmd) for example useage with parameter explanations.

Example outputs are available in [`run_example_output.zip`](examples/run_example_output.zip).

To build your own simulation setup, start from the template in [`data_sim.Rmd`](examples/data_sim.Rmd).


---

## 📊 Visualization

Example plot functions for comparing results (e.g., `plot_ssmr`, `plot_ssmd`) are provided in [`plot_ssm.Rmd`](examples/plot_ssm.Rmd).

---


## 📄 Citation

If you use this package, please cite:

> EJ Huang (2025). *SSM: A pipeline for simulating and testing phylogenetic reconstructions with state-space variation*. GitHub. https://github.com/ej91016/SSM

> (Preprint coming soon)


---

## 🪪 License

This project is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).

