# Wisdom of the crowds
This is the source code repository of the research article __"Wisdom of Crowds for Supporting the Safety Evaluation of Nanomaterials"__  by Saarimäki & Fratello et al.

## Project structure
The folder `wisdom` contains the source code used throughout the repository.

The folder `scripts` contains two scripts, one to fit the statistical response model and the other to run the cross-validation of the machine learning classifiers.

The folder `notebooks` contains the code of all the analysis performed as well as the code to reproduce the figure of the article.

The folder `data_viewer` contains the source code of the interactive data viewer.

## Data retrieval
All the data used in this project, as well as the corresponding outputs obtained by our analysis can be downloaded from this [Zenodo repository](https://doi.org/10.5281/zenodo.13884305).

To reproduce any step of the analysis, extract the dataset in this folder.

## Create the environment
To reproduce the environment used to perform the analysis, follow these instructions:

```
conda env create -n woc --file environment.yaml
conda activate woc
pip install -e .
```

## Train the models
### Response model
To fit the response model to the questionnaire data:

```python scripts/response_modeling.py```

### Machine learning classifiers
To run the cross-validation of the machine learning classifiers:

```python scripts/fit_inferred_labels.py```

This script has a number of parameters that can be changed. The default values are those used for the analysis described in the research article.

For a description of the optional parameters, refer to the script itself, or run the script with the `--help` option to print the help message.

## Data viewer
To facilitate the exploration of the questionnaire data and the analysis outputs, we developed a R shiny data viewer application. 

The easiest way to run it is to pull it from the docker hub with the following command:

```docker run --rm -it -p8080:8080 fhaive/woc_viewer```

and then open the browser to the address `localhost:8080`.

All the viewer's features are documented in the user [manual](https://github.com/fhaive/wisdom_of_the_crowds/blob/master/data_viewer/manual.pdf).
