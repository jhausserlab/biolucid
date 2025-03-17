# bioLUCID

Code for **bio**logica**L** and batch effect deconfo**U**nding by varian**C**e decompos**I**tion for scRNAseq **D**ata

## Introduction

**bioLUCID** is a Python package designed to analyze **single-cell RNA sequencing (scRNA-seq)** data by quantifying batch effects and biological variations. By applying **variance decomposition**, bioLUCID effectively separates and quantifies the sources of variation in scRNA-seq data, enabling a deeper understanding of both **technical batch effects** and **biological signals**

## Installation

```python
python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps biolucid
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
  adata = adata[adata.obs.pct_counts_mt < 5, :].copy() # the cutoff can be modified
  ```

#### Step1: Load data

- We provide an example adata containing human pbmc data sequenced on different platforms, which can be downloaded [here](https://drive.google.com/file/d/1-Uuve3sndENFDuVdSm4Ltb3Lnee7LUl8/view?usp=sharing) 

- To make sure bioLUCID works, the input data should contain:
  - Raw counts data (in `adata.X` or `adata.layers['counts']`)
  - At least 2 samples
  - In all samples, at least **2 cell types** with **more than 20 cells**
  - In all samples and in all cell types, at least **100 genes** with more than **1 UMIs** per cell

#### Step2: Initialize analyzer

- Example command is:

  ```python
  analyzer = biolucid.core.Analyzer(
      adata,
      params={
          'batch_key': 'sample',         
          'celltype_key': 'celltype',
          'min_cells': 20,
      }
  )
  ```

- To initialize analyzer, parameters are needed to pass to `Analyzer`

- All default parameters can be found in `biolucid.config.DEFAULT_PARAMS`

#### Step3: Run analyzer

- The use of bioLUCID is extremely straightforward and only requires one line of code

  ```
  analyzer.run_analysis()
  ```

- Please note: bioLUCID will call pyGSEA for enrichment calculation. The larger the database, the longer it will take. It will cost around 30 minutes for the example data.

#### Step4: Save results and visualization

- All statistical results can be summarized through `analyzer.get_key_results()` function into a dictionary, which can be further saved through `pickle`

  ```python
  results = analyzer.get_key_results()
  # import pickle
  # with open('output_file_path','wb') as f:
  #     pickle.dump(results, f)
  ```

- We also provide a data conversion function `create_analysis_df()` that can convert into a more readable dataframe (also to prepare for subsequent visualization)

  ```python
  results_df = biolucid.visualization.create_analysis_df(results)
  ```
  
- When it comes to visualization , we provided two main visualization functions: `plot_scores_bar` and `plot_scatter_analysis` 

  - `plot_scores_bar` : Creates a bar plot visualization of metrics across different batches, allowing comparison of **batch effect scores**, **normalized batch effect scores**, or **biological effect scores** with customizable display options.

    ```python
    biolucid.visualization.plot_scores_bar(results_df, score_type='Batch Effect Scores', sqrt=True, color='skyblue',figsize=(4, 6))
    ```
  - `plot_scatter_analysis` : Generates a three-panel scatter plot that simultaneously visualizes: (1) **normalized batch effects versus biological effect scores**, (2) **batch effect scores versus residual sums (with y=x reference line)**, and (3) **biological effect scores versus residual sums with Pearson correlation statistics**, providing a view of relationships between different batch effect metrics.

    ```python
    biolucid.visualization.plot_scatter_analysis(results_df)
    ```

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

- Author: Chuhanwen Sun
- Email: [silas.sun@scilifelab.se](mailto:silas.sun@scilifelab.se)
- GitHub: https://github.com/jhausserlab/biolucid