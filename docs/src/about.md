---
layout: default
title: About ccAFv2
nav_order: 2

---
# About ccAFv2
Most transcriptomic cell-cycle classifiers combine quiescent G0 and G1 cells and provide limited resolution of transitional states. We developed ccAFv2 using human neural stem cell transcriptomes to classify G1, Late G1, S, S/G2, G2/M, M/Early G1, and a quiescent-like Neural G0 state. The classifier integrates gene expression with a circular manifold representation and provides calibrated probabilities for more interpretable confidence estimates. We identified the key predictive features for each cell-cycle state and linked them to established cell-cycle programs. Benchmarking showed ccAFv2 had better performance than existing methods while resolving more states. Applications across human and mouse single-cell and single-nucleus datasets supported broad utility. Neural G0 was enriched in independently isolated quiescent neural stem cells. Analysis of embryonic cells rediscovered differentiation-associated cell-cycle regulation, demonstrating accurate classification across lineages representing all three germ layers. The higher-resolution ccAFv2 predictions were more effective at regressing out cell-cycle transcriptional signatures. Spatial transcriptomics revealed anatomically organized cell-cycle patterns in the developing mouse cortex. User-friendly ccAFv2 packages for R and Python integrate with Seurat and Scanpy workflows for single-cell and spatial studies of proliferation and neural-associated quiescence. 
![F6.large]({{ site.baseurl }}/images/F6.large.jpg)

The ccAFv2 classifier builds on the foundation established by its predecessor (ccAFv1, O'Connor et al. 2021). ccAFv2 compares well against existing cell cycle state classifiers: ccAF (v1) (O'Connor et al. 2021), Seurat (Hao et al. 2021), Tricycle (Zheng et al. 2022), SchwabeCC (Schwabe et al. 2020; Zheng et al. 2022), reCAT (Liu et al. 2017), Peco (Hsiao et al. 2020) and Cyclone (Scialdone et al. 2015). 
![F2.large]({{ site.baseurl }}/images/F2.large.png)

To improve machine learning classification of cell cycle states, we developed a probabilistic approach to represent each cell�s position along the circular cell cycle manifold. By encoding angular positions as probability distributions across equally sized bins, this method preserves circular continuity while accounting for uncertainty in cell positioning. We evaluated multiple binning resolutions and found that 10 angular bins provided an optimal balance of resolution, model stability, and interpretability. Increasing the resolution to 36 bins did not significantly improve classification performance, supporting the selection of 10 bins for the final ccAFv2 classifier.
![F3.large]({{ site.baseurl }}/images/F3.large.png)

The core algorithm of the ccAFv2 is a fully connected artificial neural network (ANN) implemented using the PyTorch (v2.2.1) to construct ANNs. It takes the expression levels of 861 highly variable genes and 10 angular bins as input and classifies cells into seven distinct cell cycle states: Neural G0, G1, Late G1, S, S/G2, G2/M, and M/Early G1. The model incorporates dropout layers to prevent overfitting and uses Stochastic Gradient Descent (SGD) for optimization. The result is a powerful, flexible tool for analyzing cell cycle dynamics.
![F1.large]({{ site.baseurl }}/images/F1.large.png)

### References 
O'Connor SA, Feldman HM, Arora S, Hoellerbauer P, Toledo CM, Corrin P, Carter L, Kufeld M, Bolouri H, Basom R, et al. 2021. *Neural G0: a quiescent-like state found in neuroepithelial-derived cells and glioma.* Mol Syst Biol 17: e9522. [PubMed](https://pubmed.ncbi.nlm.nih.gov/34101353/)

Hao Y, Hao S, Andersen-Nissen E, Mauck WM, Zheng S, Butler A, Lee MJ, Wilk AJ, Darby C, Zager M, et al. 2021. *Integrated analysis of multimodal single-cell data.* Cell 184: 3573–3587.e2 [PubMed](https://pubmed.ncbi.nlm.nih.gov/34062119/)

Schwabe D, Formichetti S, Junker JP, Falcke M, Rajewsky N. 2020. *The transcriptome dynamics of single cells during the cell cycle.* Mol Syst Biol 16: e9946.[PubMed](https://pubmed.ncbi.nlm.nih.gov/33205894/)

Zheng SC, Stein-O’Brien G, Augustin JJ, Slosberg J, Carosso GA, Winer B, Shin G, Bjornsson HT, Goff LA, Hansen KD. 2022. *Universal prediction of cell-cycle position using transfer learning.* Genome Biol 23: 41. [PubMed](https://pubmed.ncbi.nlm.nih.gov/35101061/)

Liu Z, Lou H, Xie K, Wang H, Chen N, Aparicio OM, Zhang MQ, Jiang R, Chen T. 2017. *Reconstructing cell cycle pseudo time-series via single-cell transcriptome data.* Nat Commun 8: 22. [PubMed](https://pubmed.ncbi.nlm.nih.gov/28630425/)

Hsiao CJ, Tung P, Blischak JD, Burnett JE, Barr KA, Dey KK, Stephens M, Gilad Y. 2020. *Characterizing and inferring quantitative cell cycle phase in single-cell RNA-seq data analysis.* Genome Res 30: 611–621. [Full Text](https://genome.cshlp.org/content/30/4/611?ijkey=b6e258cf40d7384bb32fc892319f1bee12f61ae9&keytype2=tf_ipsecsha)

Scialdone A, Natarajan KN, Saraiva LR, Proserpio V, Teichmann SA, Stegle O, Marioni JC, Buettner F. 2015. *Computational assignment of cell-cycle stage from single-cell transcriptome data.* Methods 85: 54–61 [PubMed](https://pubmed.ncbi.nlm.nih.gov/26142758/)
