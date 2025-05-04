# inferCNAsc

**Version: 0.1 (development)**  
A Python package for inferring copy number alterations (CNAs) from single-cell RNA-sequencing (scRNA-seq) data.

---

## Overview

**inferCNAsc** enables inference of DNA copy number gains and losses directly from single-cell gene expression data. Designed for use with `AnnData` objects, it supports applications in stem cell biology, cancer research, and genome stability analysis. It is particularly optimized for identifying CNAs in human pluripotent stem cells (PSCs) and their derivatives.

---

## Features

- Inference of CNAs from scRNA-seq count data
- Cell-level or group-level CNA assignments
- Visualization tools for chromosome-wide CNA footprints
- Benchmarking utilities for simulated data validation
- Compatible with `scanpy`, `anndata`, `pandas`, and other ecosystem tools

---

## Installation

This project is under active development. Clone the repository and install in editable mode:

```bash
git clone https://github.com/YOUR_USERNAME/inferCNAsc.git
cd inferCNAsc
pip install -e .
```

## Authors

This project was co-developed, in equal collaboration, by:

- **Alejandro J. Soto Franco**
- **Raeann Kalinowski**
- **Amy Liu**

Computational Stem Cell Biology Final Project, Spring 2025

> Johns Hopkins University Department of Biomedical Engineering

---

## Citation

If you use `inferCNAsc` in your work, please cite:

> Soto Franco A.J., Kalinowski R., Liu A. *inferCNAsc: a Python toolkit for copy number inference from single-cell transcriptomes*. In preparation (2025)
