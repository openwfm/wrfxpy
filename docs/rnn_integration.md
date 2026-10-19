# Integrating Recurrent Neural Networks for Fuel Moisture Prediction into FMDA

## Clone the RNN Compatible Branch 

```bash
git clone https://github.com/openwfm/wrfxpy
cd wrfxpy
git switch dev_fmda_JH
```

## Get static data
Acquire the elevation static data from here:

```bash
mkdir static
cd static
wget https://demo.openwfm.org/web/wrfx/ds.terrainh.bin
tar xvfj hrrr_terrain.tbz
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

The import check below also requires `xgboost` and `pytables`:
```bash
conda install -c conda-forge xgboost pytables
```
Use the scikit-learn version recorded in the pretrained scaler. The supplied
operational scaler used in the October 19, 2026 validation was saved with
scikit-learn 1.7.1.

Then, navigate to the project root for `ml_fmda` and build the local package
```bash
cd ml_fmda
pip install -e .
```

Confirm instillation at this point 
```bash
python -c "import numpy, pandas, sklearn,  tensorflow; print('numpy', numpy.__version__); print('tensorflow', tensorflow.__version__)"
```


### 3. Confirm the Package Installation

From the `ml_fmda` project root, run its small training and prediction script:

```bash
cd /path/to/ml_fmda
python tests/test_rnn.py
```

At `ml_fmda` revision `89220ac`, this file places its assertions inside
`if __name__ == '__main__'`. Invoking pytest on it collects no tests and exits
with code 5. Run the wrfxpy tests after configuring the model assets below.


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

### Test the wrfxpy Integration

`src/tests/test_rnn.py` reads `etc/rnn_cycler.json`, model parameters, the scaler,
weights, and transfer table during test collection. These inputs must exist
before invoking pytest, even though package installation does not need them.

From the wrfxpy project root:
```bash
python -m pytest -v src/tests/test_rnn.py
python -m pytest -v src/tests/test_rnn.py src/tests/test_rnn_integration.py src/tests/test_rnn_cycler.py
```

The additional tests check nonsquare grids, WPS readback, blending, feature
scaling, prediction batching, and synthetic cycler continuation. They also
check that incomplete output is not treated as a completed cycle. They require
`etc/conf.json`; an empty JSON object `{}` is sufficient in an isolated test
checkout. Keep an operational checkout's existing configuration.

The synthetic cycler test supplies weather arrays, returns no observations,
and disables rendering. It does not validate HRRR ingestion, assimilation,
forecast accuracy, or rendered products. See `notes/rnn-validation-20261019.txt`
for the observed failures and remaining real-data checks.
 
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

