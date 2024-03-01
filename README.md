MS CITE-seq analysis
---

* [Data integration](step1_merge_data_qc.html)

* [Cell type identification](step2_celltype_annotation.html)

* [DEG at the top level](step3_deg_analysis.html)

* [Subtype identification](step4_subtype_analysis.html)

* [DEG at the subtype level](step5_deg_subtype_analysis.html)

# Methods

## Batch correction and visualization of scRNA-seq data 

Our single-cell data set consists of five different batches of CITE-seq experiments [PMID: 28759029]. Each CITE-seq experiment was also multiplexed to capture general CD4+ T cells and CD25+ enriched CD4+ T cells. In total, we profiled expression vectors of 33,578 features (surface protein antibodies and genes) across 41,957 cells; overall, the data matrix was sparse, containing only 4% of non-zero elements (n = 57,290,814), implicating many features potentially dropped out due to technical limitations.

### Quality Control (QC) of cells

We identified cells that unexpectedly contain a low number of non-zero gene expressions and highly affected mitochondrial mRNA. Lowly-expressed cells often implicate technical issues in droplet-based single-cell sequencing technology, whereas high mitochondrial concentrations can be a good sign of dying, or burst, cells and droplets. Instead of borrowing threshold values from previous unrelated studies, we first determined the cutoffs for the highest possible mitochondrial fractions among viable cells by performing k-means clustering implemented in R with three components and removed cells that fell in the highest group; we also identified cells that expressed too few gene expressions by k-means clustering with the two components and removed cells in the low-expressed cluster. We visually confirmed bimodal distributions in both types of statistics, and cells retained belong to major groups not to eliminate too many cells. We retained 35,488 cells and discarded 6,462 cells. We further tested whether the retained cells were well distributed within the principal component (PC) space. We used the top five PCs computed on the log1p-normalized gene expression matrix and identified a small set of cells exhibiting a bimodal distribution along the first and second PC axes. We conducted two-component k-means clustering to identify and remove one outlier group of cells. 

### Quality Control (QC) of features

Lowly-expressed genes can induce outlier cell groups and may introduce unwanted biases and scattered clustering patterns. We computed each gene's statistics, including the number of cells expressed, the mean and standard deviation (SD), and the coefficient of variation (CV; the mean value divided by the SD). We largely identified two distinctive clusters in a two-dimensional density plot of the CV and the number of non-zero cells for each gene. We removed a gene if it is only expressed in less than 500 cells or has a CV value of less than 1.25 (not much variation across cells).

## Cell type assignment

Overall, we performed cell-type annotation in two stages. First, we split cells into four major cell types--mTreg, nTreg, mTconv, and nTconv. Later, among the mTreg and mTconv cell groups, independently, we resolved cell clusters in unsupervised graph-based clustering methods (Leiden) and called cell clusters that correspond to subtypes within each major cell type.

### Marker-based annotation

In the first stage, we first established cell groups based on known protein surface marker activities. We systematically sorted cells according to the following decision rule, which was previously suggested in cell sorting experiments [PMID: 30449622]:

* mTreg: CD3+, CD4+, CD8-, CD25+/CD127-, CD45RA-/CD45RO+
* mTconv: CD3+, CD4+, CD8-, CD25-/CD127+, CD45RA-/CD45RO+
* nTreg: CD3+, CD4+, CD8-, CD25+/CD127-, CD45RA+/CD45RO-
* nTconv: CD3+, CD4+, CD8-, CD25-/CD127+, CD45RA+/CD45RO-

We define a latent indicator variable $z_{jk}$ to mark the assignment of a cell $j$ to a cell type $k$ and estimate the posterior probability of $z_{jk}=1$ by the stochastic expectation maximization (EM) algorithm. We assume that a normalized vector for each cell $\mathbf{x}_{j}$ follows von Mises-Fisher (vMF) distribution [https://doi.org/10.1098/rspa.1953.0064] with cell type k-specific mean vector $\mu_{k}$ and shared concentration parameter $\kappa$: $P(\mathbf{x}_{j}|z_{jk}=1) \sim \exp⁡(\kappa \theta^{\top} \mathbf{x}_{j})$. We modified the existing EM algorithm [PMID: 15384557] and enforced the sparsity of the mean vector $\mu$ using the prior knowledge of the cell-type-specific activity of marker proteins/genes/features. Simply, we allow $\mu_{gk}$ to take a non-zero value if and only if a feature $g$ is a known marker for the cell type $k$. This stochastic EM algorithm is implemented in `rcpp_mmutil_annotate_columns` function made publicly available within our `mmutilR` R package. 

### Refinement of the initial cell type assignment

Additionally, we performed graph-based clustering to visualize and confirm the estimated cell types. After removing scattered singletone clusters, we were able to partition cells into nine distinctive clusters. To avoid putative batch-specific clusters, we performed batch-balancing k-nearest neighbour (BBKNN) analysis [CITE] and kept cell pairs/edges only if both endpoints/cells mutually included the other cell among k-nearest neighbour cells. We found nearest neighbours by a scalable approximate neighbourhood search algorithm, and the distances between cells were calculated by the Euclidean distance of the top 50 PCs.

Although gene expression and surface marker antibody activities generally agree with each other, we regretfully spotted sporadic discrepancies between them. Therefore, we investigated whether a cell cluster contains cells of the same cell type based on the surface marker. If a cell cluster consists of mixed cell types, we eliminated putative doublet-like clusters to ensure that the cell type annotation results are intact not only globally but also locally. We retained a cell cluster if more than half of the cells were annotated to the same cell type.

### Subtype analysis 

We resolved cell clusters within mTreg and mTconv populations independently and identified distinctive subtypes within each major T cell population. Since raw sequencing profiles are not robust enough, we projected cells sampled from the original transcriptomic space of 12k+ features into a thirty-dimensional latent space. We conducted randomized singular value decomposition (SVD) on log1p-transformed data. We did not find any significant batch-specific patterns for mTconv cells. However, we found an additional group of outliers among mTreg cells, specifically deviating in the PC2 axis, so these cells were removed, and fifty PCs were reestimated using randomized SVD. We then controlled and removed putative batch-specific effects by constructing BBKNN graphs, searching with 50 nearest neighbours for each cell across five different batches in the PC space. We ran a graph-based Leiden clustering method on the BBKNN graphs and identified three clusters in both mTconv and mTreg cells.

