# bioLUCID

Code for bioLUCID: quantifying batch and biological effect in single-cell transcriptomics by hypothesis-driven variance decomposition

## Introduction

bioLUCID is a Python package specifically designed for single-cell RNA sequencing (scRNA-seq) data analysis, focused on quantification of batch effects and biological variations. Through variance decomposition, bioLUCID effectively distinguishes and quantifies variations from different sources (samples), to better understand technical batch effects and biological signals within scRNA-seq data.

## Installation

1. **System Requirements**

- Python >= 3.8

2. **pip installation**

```python
pip install -i https://test.pypi.org/simple/ biolucid==1.0.6
```

## Quick Start

- See `tests/Human_PBMC_data_test.ipynb` to walk through bioLUCID tutorial.

#### Step0: Preprocessing

- In the tutorial with provided example data, the preprocessing step has been done, including quality control and cell type annotation. 

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

#### Step1: Load data

- We provide an example adata containing human pbmc data sequenced on different platforms, which can be downloaded [here](https://drive.google.com/file/d/1-Uuve3sndENFDuVdSm4Ltb3Lnee7LUl8/view?usp=sharing)
- To make sure bioLUCID works, the input data should contain:
  - Raw counts data (in `adata.X` or `adata.layers['counts']`)
  - At least 2 samples
  - In all samples, at least **2 cell types** with **more than 20 cells**
  - In all samples and in all cell types, at least **100 ubiquitously expressed genes** with more than **1 UMIs** per cell

#### Step2: Initialize analyzer

- Example command is:

  ```python
  analyzer = biolucid.core.BatchEffectAnalyzer(
      adata,
      params={
          'batch_key': 'sample',         
          'celltype_key': 'celltype',
          'min_cells': 20,
      }
  )
  ```

- To initialize analyzer, parameters are needed to pass to `BatchEffectAnalyzer`

- All default parameters can be found in `biolucid.config.DEFAULT_PARAMS`

#### Step3: Run analyzer

- The use of bioLUCID is extremely straightforward and only requires one line of code

  ```
  analyzer.run_analysis()
  ```

#### Step4: Save results and visualization

- All statistical results can be summarized through `analyzer.results` function into a dictionary, which can be further saved through `pickle`

  ```python
  results = analyzer.results
  # import pickle
  # with open('output_file_path','wb') as f:
  #     pickle.dump(results, f)
  ```

- We also provide a data conversion function `results_to_df()` that can convert into a more readable dataframe (also to prepare for subsequent visualization)

  ```python
  global_results_df,per_sample_results_df = biolucid.visualization.results_to_df(results)
  ```

- When it comes to visualization , we provided one visualization function: `plot_scatter_analysis` to generate qsh_qsp plot

    ```python
    biolucid.visualization.plot_scatter_analysis(per_sample_results_df)
    ```

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

- Author: Chuhanwen Sun
- Email: [silas.sun@scilifelab.se](mailto:silas.sun@scilifelab.se)