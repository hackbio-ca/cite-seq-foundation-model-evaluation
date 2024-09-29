# Surface protein abundance prediction from RNA expression - an analysis and benchmark of foundation models

 <img src="static/citeseq.png" height="250px" align="right"/>

#### Toronto Bioinformatics Hackathon 2024 - Philip Fradkin, Cait Harrigan, Hassaan Maan

We aim to evaluate the mutual information between RNA expression and surface protein abundance through a large CITE-seq dataset of peripheral blood mononuclear cells (PBMCs). CITE-seq combines these two measurements such that RNA expression and surface protein abundance measurements are sampled for the same cells. Further, we evaluate foundation models on CITE-seq data to determine if zero-shot representations in these models can allow for more accurate prediction of protein abundance from RNA expression. Our analysis offers a unique perspective on this prediction task by determining real bounds on mutual information and evaluating the performance of cutting-edge foundation models on this task. 

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Abstract

The relationships between the different levels of information in the central dogma hold important answers in our understanding of biology and disease. Although multi-modal measurements for sequencing data have become increasingly common, the mutual information and covariation between these modalities is often difficult to quantify and varies based on cell and tissue context. Previous approaches to predicting protein expression from RNA have relied on specialized models that transform the space between one modality to another. However, this does not answer fundamental questions about their relationships. In this work, we aim to explicitly determine information bounds between RNA and surface protein expression through a large CITE-seq dataset of peripheral blood mononuclear cells (PBMCs) that jointly measures both modalities in single-cells. After determining that mutual information can be improved through incorporating foundation model embeddings, we evaluate both single-cell and sequence-based foundation models in the surface protein prediction task. Our analysis offers a unique perspective on this prediction task by determining real bounds on mutual information and evaluating the performance of cutting-edge foundation models in a novel scenario with real-world applications, such as the study of disease mechanisms and developing RNA-based therapeutics and diagnostics.

## Installation

Conda is necessary to install the necessary dependencies for this project. To install Conda, follow the instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

To install the necessary libaries and dependencies, run the following commands after installing conda:

```bash
# Install the dependencies from the env.yaml file
conda env create -f env.yaml

# Activate the environment
conda activate cite_seq
```

## Quick Start

As the CITE-seq data needs a lot of processing, and the embeddings for scGPT, Orthrus and ESM2 are ardous to generate, we provide a MuonData file that contains the gene expression (SCT), surface protein abundance (ADT), scGPT embeddings, Orthrus embeddings, and ESM2 embeddings. The original CITE-seq data can be found [here](https://atlas.fredhutch.org/nygc/multimodal-pbmc/) in the downloads section.

The fully processed MuonData object can be downloaded from [here](https://drive.google.com/file/d/1H8L7eaLkcHthNa7McQKfAQLrHwkFxOIO/view?usp=sharing).

To run the analysis notebooks for mutual information, notebooks 04-06 can be run in order. A GPU is recommended and the `cite_seq` environment should be used. 

```python
# For example
jupyter notebook notebooks/04_mi_test_citeseq_adt_sct.ipynb
```

Other analysis notebooks in `notebooks` can be run in a similar manner.

To run the experiments for the foundation models for prediction of surface protein (ADT) from gene expression (SCT), a configuration file from `configs` can be used:

```python
python train.py --config_file [config_file.yaml] \
    --seed 3023 \
    --linear_model=true \
    --note_suffix="_linear" \
```

We recommend using a GPU for training the models. The configuration files can be found in the `configs` folder.

## Contribute

Contributions are welcome! If you'd like to contribute, please open an issue or submit a pull request. See the [contribution guidelines](CONTRIBUTING.md) for more information.

## Support

If you have any issues or need help, please open an [issue](https://github.com/hackbio-ca/cite-seq-foundation-model-evaluation/issues) or contact the project maintainers.

## License

This project is licensed under the [MIT License](LICENSE).
