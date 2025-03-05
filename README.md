# Code for the analysis in the paper `Recent gains in global terrestrial carbon stocks are mostly stored in nonliving pools`

All data to support the findings of the study are publicly available and could be downloaded using the `00a_download_data.ipynb` notebook. The only data not publicly available is the L-VOD biomass data, which could be received through a request to the authors.

The code was tested on python 3.11.9

To run the code, we recommend installing [uv](https://github.com/astral-sh/uv), setting up a new python environment and installing the required dependencies using the `requirements.txt` file:


```bash
uv venv --python 3.11.9
source .venv/bin/activate
uv pip install -r requirements.txt
```

You will also need to install GDAL using conda. First, make sure conda or miniconda is installed on your system. Then run the following command:

```bash
conda install -c conda-forge gdal
```

To run the code, first download the data using the `00a_download_data.ipynb` notebook. The only data not publicly available is the L-VOD biomass data, which could be received through a request to the authors.

Once all of the data has been downloaded, you can run the scripts in consecutive order to produce the results and figures. The scripts are numbered in the order they should be run.