# bioLUCID

bioLUCID is a Python package to quantify batch and biological effects in a collection of single-cell RNA sequencing (scRNAseq) samples.

Batch-correcting single-cell gene expression data enables multi-sample analyses but risks erasing sample-specific biology along with batch effects. To address this, we introduce bioLUCID, a heuristic to decompose gene expression variation into biological and batch effects. This decomposition allows (i) integrating cell embeddings across samples without erasing sample-specific biology, and (ii) flagging samples with excessive batch effects, with benefits for multiple applications such as exploring transcriptional heterogeneity and multi-omics integration.

A pre-print describing the approach is available on [Research Square](https://www.researchsquare.com/article/rs-7290456/v1)

Silas Sun, Alper Eroglu, Jean Hausser at Karolinska Institutet & SciLifeLab, Stockholm Sweden

## Quick start

```python
$ pip install -i https://test.pypi.org/simple/ biolucid==1.1.0
```

```python
# Replace the line below with your adata object
adata = ...

# Run bioLUCID
biolyzer = biolucid.core.BatchEffectAnalyzer(
    adata,
    params={'batch_key': 'sample',
            'celltype_key': 'celltype'})
biolyzer.run_analysis()

# Visualize the results
global_results_df, per_sample_results_df =
biolucid.visualization.results_to_df(biolyzer.results)
[per_sample_table]
biolucid.visualization.plot_scatter_analysis(per_sample_results_df)
[plot] 
```


# User manual

## Installation

1. **System Requirements**

- Python >= 3.10

2. **pip installation**

```python
pip install -i https://test.pypi.org/simple/ biolucid==1.1.0
```

## Quick Start

- See `tests/Human_PBMC_data_test.ipynb` to walk through bioLUCID tutorial.

### Step0: Preprocessing

- bioLUCID requires data in AnnData format (adata object), which is the standard format for single-cell analysis in Python.

- In the tutorial with provided example adata, the preprocessing step has been done, including quality control and cell type annotation. 

- However, if the data is raw, we recommend a standard QC step before running bioLUCID:

  ```python
  # filter cells and genes
  sc.pp.filter_cells(adata, min_genes=200)
  sc.pp.filter_genes(adata, min_cells=3)
  
  # filter MT genes
  adata.var['mt'] = adata.var_names.str.startswith('MT-')
  sc.pp.calculate_qc_metrics(
      adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True
  )
  adata = adata[adata.obs.pct_counts_mt < 20, :].copy() # the cutoff can be modified
  ```

### Step1: Load data

- The general guideline is to give to bioLUCID the set of samples that you plan to contrast in the downstream analysis. Avoid processing sets of samples from very different tissues tissue types. For example, passing whole tissue scRNAseq from lung and brain to bioLUCID will cause practical issues (finding common cell types and ubiquitously-expressed genes) as well as fundemental ones (all lung cell types may have common transcriptional differences compared to all brain cell types, which bioLUCID will interpret as batch).
- We provide an example adata containing human pbmc data sequenced on different platforms, which can be downloaded [here](https://drive.google.com/file/d/1-Uuve3sndENFDuVdSm4Ltb3Lnee7LUl8/view?usp=sharing)
- To make sure bioLUCID works, the input data should contain:
  - unnormalized UMI counts data in `adata.X` or `adata.layers['counts']`
  - sample identifiers in `adata.obs`, stored as a column `batch` by default (can be parametrized, see below)
  - cell type strings in `adata.obs`, stored as a column `celltype` by default (can be parametrized, see below)

These are the minimal data requirements to run bioLUCID:
- At least 2 samples
- In all samples, at least **2 cell types** with **more than 20 cells each**
- In all samples and in all cell types, at least **100 ubiquitously expressed genes** with more than **1 UMIs** per cell on average

### Step2: Initialize analyzer

- Example command is:

  ```python
  analyzer = biolucid.core.BatchEffectAnalyzer(
      adata,
      params={
          'batch_key': 'sample',         
          'celltype_key': 'celltype',
          'min_cells': 20,
          'abundant_gene_threshold': 1,
          'min_abundant_genes': 100,
      }
  )
  ```

- To initialize analyzer, parameters are needed to pass to `BatchEffectAnalyzer`

- All default parameters can be found in `biolucid.config.DEFAULT_PARAMS`

### Step3: Run analyzer

- The use of bioLUCID is extremely straightforward and only requires one line of code

  ```python
  analyzer.run_analysis()
  ```

### Step4: Save results and visualization

- All statistical results can be summarized through `analyzer.results` function into a dictionary, which can be further saved through `pickle`

  ```python
  results = analyzer.results
  # import pickle
  # with open('output_file_path','wb') as f:
  #     pickle.dump(results, f)
  ```

- We also provide a data conversion function `results_to_df()` that can convert into a more readable dataframe (also to prepare for subsequent visualization)

  ```python
  global_results_df, per_sample_results_df = biolucid.visualization.results_to_df(results)
  ```

#### Results Example

The `per_sample_results_df` will contain results in the following format:

|  | q_sh_score_per_batch | q_sp_score_per_batch | b_score_per_batch | recommendation |
|--------|---------------------|---------------------|------------------|----------------|
| 10x-Chromium-v2-A | 0.160767 | 0.075598 | 0.818920 | Drop |
| 10x-Chromium-v2-B | 0.173909 | 0.054345 | 0.911036 | Drop |
| 10x-Chromium-v3 | 0.240570 | 0.101513 | 0.848854 | Drop |

**Column Descriptions:**
- `q_sh_score_per_batch`: Estimated magnitude of variation shared across the cell types
- `q_sp_score_per_batch`: Estimated magnitude of variation specific to a cell type
- `b_score_per_batch`: The relative contribution of batch effects to gene expression
- `recommendation`: Recommendation for sample selection

Both q_sh and q_sp are measured in units of standard deviation on log expression. See our manuscript for details. At the same time, we also provide advice on the selection of each sample.

#### Visualization

When it comes to visualization, we provided one visualization function: `plot_scatter_analysis` to generate qsh_qsp plot

```python
biolucid.visualization.plot_scatter_analysis(per_sample_results_df)
```

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

- Author: Chuhanwen Sun
- Email: [silas.sun@scilifelab.se](mailto:silas.sun@scilifelab.se)
