# Integrating Recurrent Neural Networks for Fuel Moisture Prediction into FMDA

## Clone the RNN Compatible Branch 

```python
git clone https://github.com/openwfm/wrfxpy
cd wrfxpy
git switch dev_fmda_JH
```

## Install Necessary Packages

The RNN functionality is built in `tensorflow`, with custom model classes in a git project `ml_fmda`. The `ml_fmda` package is installed locally after cloning. 

### 1. Installing the `ml_fmda` package

The `ml_fmda` git project is a larger code base with support for data retrieval and model training. Inside this project is the package with core functionality: `ml_fmda/src/ml_fmda`. For use within `wrfxpy`, first clone the main project repo

```bash
git clone https://github.com/openwfm/ml_fmda.git
cd ml_fmda
git checkout main 
```


### 2. Build the Conda Environment

It is strongly recommended you build this environment from scratch, rather than appending `tensorflow` to an existing conda environment that works with `wrfxpy`. 

TODO: test a more elegant build recipe, but this works as of Oct 14 2026 so I don't want to touch it

```bash
conda create -n wrfx_ml python=3.12
conda activate wrfx_ml

conda install -c conda-forge gdal
conda install -c conda-forge netcdf4 h5py pyhdf pygrib simplekml
conda install -c conda-forge f90nml lxml pytz pandas scipy
conda install -c conda-forge basemap paramiko dill psutil flask wgrib2
conda install -c conda-forge toml synopticpy==2024.9.0
conda install -c conda-forge scikit-learn pyyaml joblib tensorflow
pip install python-cmr shapely
conda install -c conda-forge pyarrow
conda install pytest
```

Then, navigate to the project root for `ml_fmda` and build the local package
```bash
cd ml_fmda
pip install -e .
```

Confirm instillation at this point 
```bash
python -c "import numpy, pandas, sklearn, xgboost, tables, tensorflow; print('numpy', numpy.__version__); print('tensorflow', tensorflow.__version__)"
```


### 3. Confirm Installation with PyTest

To confirm the `ml_fmda` install, from the `ml_fmda` project root run the tests 

```bash
cd /path/to/ml_fmda
pytest tests/test_rnn.py
```

To confirm the integration into `wrfxpy`, from the `wrfxpy` project root run the tests
```bash
cd /path/to/wrfxpy
pytest src/tests/test_rnn.py
```


## Get the Pretrained RNN Model Weights

Integration into `wrfxpy` assumed pretrained RNN model weights. A set of pretrained model weights and associated scalers/configuration are distributed separately from the source package as a ZIP archive: `fm_transfer.zip`. The archive includes a README and the models from the training and transfer-learning workflows. Two main directories store the model weights, parameter files for building the model architecture, and information about how they were trained. 

```text
fm_transfer/
|-- README.md
|-- fm_target_models/
|-- fm10_source_model/
```

After obtaining the archive, extract it to a local directory:
```bash
unzip fm_transfer.zip
```

Once the pretrained model weights are extracted, update a configuration file that tells `wrfxpy` where to find model weights
```bash
cp etc/rnn_cycler.json.initial etc/rnn_cycler.json
```

Update the paths
```bash
    "rnn_model_path": "/path/to/fm10_source_model",
    "transfer_path": "/path/to/fm_target_models"
```
 
## Run FMDA with RNN integration

Analysis mode:
```bash
./hrrr_cycler.sh a
```

Forecast mode:
```bash
./hrrr_cycler.sh f
```

If successful, model outputs should be found within the parent output directory `fmda-DOMAIN-YYYYMMDD-HH`, with the following naming conventions:
* RNN 
	- Netcdf format: `rnn-DOMAIN-YYYYMMDD-HH.nc` 
	- Geogrid format: `rnn-DOMAIN-YYYYMMDD-HH.geo`
	- WPS format: `RNN:DOMAIN-YYYYMMDD-HHfHH`
* Model Blend  
        - Netcdf format: `blend-DOMAIN-YYYYMMDD-HH.nc`                         
        - Geogrid format: `blend-DOMAIN-YYYYMMDD-HH.geo`
        - WPS format: `BLEND:DOMAIN-YYYYMMDD-HHfHH`

Postprocessing outputs should be found within the parent directorry `YYYYMM/fmda-DOMAIN-YYYYMMDD-HH`, with the naming conventions:
* Time lag model: `1-hr DFM-raster.png`
* RNN:            `1-hr DFM_RNN-raster.png`
* Model Blend     `1-hr DFM_BLEND-raster.png`

## AI use disclosure
OpenAI Codex and ChatGPT were used as coding assistants and documentation editors. The author reviewed and directed their use and remains responsible for the software and documentation.

