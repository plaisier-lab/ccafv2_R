---
title: Main
layout: home
nav_order: 1
---
# ccAFv2 R package 
![Logo]({{ site.baseurl }}/images/logo.png)

**ccAFv2** is a high-resolution classifier for cell cycle state identification in single-cell RNA-seq (scRNA-seq) data. It was trained on human neural stem cells and classifies six distinct cell cycle states (G1, Late G1, S, S/G2, G2/M, and M/Early G1) along with a Neural G0 state (G0). ccAFv2 incorporates a tunable parameter to refine classification certainty and outperforms existing methods while supporting additional cell cycle states.

Key Applications of ccAFv2:
- Identify cell cycle states in diverse cell types and germ layers.
- Classify cells, nuclei, and spatial transcriptomics data across species.
- Regress cell cycle expression patterns to reveal biological signals.
- Use as an [R package](https://github.com/plaisier-lab/ccafv2_R) (Seurat integration) or [PyPI package](https://pypi.org/project/ccAF/) (scanpy integration).

# Citing ccAFv2
If you include or rely on ccAF2 when publishing research, please adhere to the following citation guide:

**Citation for ccAFv2 (version 2):**

*Classifying cell cycle states and a quiescent-like G0 state using single-cell transcriptomics.* Samantha A. O‘Connor, Leonor Garcia, Rori Hoover, Anoop P. Patel, Benjamin B. Bartelle, Jean-Philippe Hugnot, Patrick J. Paddison, Christopher L. Plaisier. bioRxiv [Preprint]. 2024 Apr 20:2024.04.16.589816. doi: 10.1101/2024.04.16.589816. PMID: 38659838

**Citation for ccAF (version 1):**

*Neural G0: a quiescent-like state found in neuroepithelial-derived cells and glioma.* Samantha A. O'Connor, Heather M. Feldman, Chad M. Toledo, Sonali Arora, Pia Hoellerbauer, Philip Corrin, Lucas Carter, Megan Kufeld, Hamid Bolouri, Ryan Basom, Jeffrey Delrow, Jose L. McFaline-Figueroa, Cole Trapnell, Steven M. Pollard, Anoop Patel, Patrick J. Paddison, Christopher L. Plaisier. Mol Syst Biol. 2021 Jun;17(6):e9522. doi: 10.15252/msb.20209522. PMID: 34101353

{: .Acknowledgments-title }
>**Acknowledgments**
>
>Samantha A. O‘Connor developer since 2021, Christopher L. Plaisier developer since 2021, maintainer 
>
>This website was developed by Rori Hoover and the Plaisier Lab using [Just the Docs](https://github.com/just-the-docs/just-the-docs).
>
> For other great packages from the Plaisier Lab, please check [here](https://github.com/plaisier-lab).

# Support 
Found a bug or would like to see a feature implemented? Feel free to submit an [issue](https://github.com/plaisier-lab/ccAFv2_R/issues/new). Your help to improve ccAF is highly appreciated. 



