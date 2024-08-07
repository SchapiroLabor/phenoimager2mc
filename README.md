# phenoimager2mc
Formatting PhenoImager .tif output files to be compatible with the MCMCIRO pipeline and ASHLAR.

## Description
**Raw data:**
The PhenoImager software outputs one float32 .tif file per tile and cycle containing all channels. The metadata is unstandardized. 
**Goal:**
For analysing the data within MCMCIRO (initial registration with ASHLAR), one stacked ome-tif file per channel containing all tiles and cycles is required. 
**Steps in this module:**
* Extraction of metadata from unstandardized tif files
* Creation of stacked and correct ome-tiff files readable for ASHLAR
* Conversion from float32 to uint16
* Normalization to max or 99th percentile (user's choice)

## Usage

### CLI
#### Input
The CLI script `scripts/phenoimager2mc.py` requires 3 inputs
* The path to the folder containing all .tif files from one cycle with `-i` or `--input`
* The number of markers that was used in this cycle with `-m` or `--num_markers`
* The normalization method that the intensities per cycle should be normalized with. Either 99th or max with `-n` or `--normalization`

#### Output
* Output .tif file containing all tiles and channels of one cycle with `-o` or `--output`

### Docker usage

If you want to run the module directly from a pre-configured container with all the required packages, you can either build the docker container yourself or pull it from the Github container registry.

To build the container run:

```
git clone https://github.com/SchapiroLabor/phenoimager2mc.git
docker build -t phenoimager2mc:latest .
docker run phenoimager2mc:latest python phenoimager2mc.py
```

To pull the container from the Github container registry (ghcr.io):

```
## Login to ghcr.io
docker login ghcr.io

## Pull container
docker pull ghcr.io/schapirolabor/phenoimager2mc:latest
```
