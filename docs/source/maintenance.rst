Maintenance
===========
This page is about how to generate dataset, produce new library.

Dataset Info
------------

.. note::
    This part is adapted from `document <https://docs.google.com/document/d/1SGd7ynyIR5mS4uTZtaV2XU5f1GMde6DjMPQBccf4bv8/edit?usp=sharing>`_ written by Esther

Each recording data “group” is governed by one excel document. Currently we are using `chamLabRecordingVariables.xlsx <https://drive.google.com/drive/folders/1fpUvxR17hc5CaAnXwgyjzDOEguGLr4Bh?usp=sharing>`_

Each excel document has a sheet titled “datasetInfo”, this sheet is responsible for:

* Defining all necessary variables for all datasets used within the data “group”
* Of the variables defined, the following must be included 1: “fileName”, 2: “RAWfilePath”, 3: “samplingRate”, 4: “dataChannels”, 5: “cacheDirPath”, 6: “recordingLength” (this one may not be necessary in the future), with the naming as written
* The data type of each variable entry for a single dataset needs to be specified: integers = “int”, floats = “float”, strings = “str”, a list of integers = “intList”, list of floats = “floatList”, list of strings = “strList”

Each excel document also has a sheet titled “dataCacheFileNames”, this sheet is responsible for titling the files containing data arrays/structure that will be cached, such as filtered raw data, spike-sorting outputs, spiking times, pca outputs, etc. These file names will be the naming convention for all dataset within the data “group”:

* Naming convention here should allow for specification of data analysis parameters.  It is important for files containing data to have all important information about the configuration and parameters that yield the data in the filename.  In the case that the file is copied outside of the cache, important parameters should be listed in the title.
* Of the variables defined, the following must be included 1: “pca_scores”, 2: “amplitudes”, 3: “waveforms” 4: “spikes2”, 5: “spikes we should make better naming convention for spikes and spikes2 6: “recording_cmr_filename”, 8: “mountainSortExport”, 9: ”exportMatlabData”

Each excel document also has a sheet titled “imageCacheFileNames”, this sheet is responsible for titling the files containing image outputs of the data analysis, such as unit waveforms, histograms, etc. All files will be saved as .png files:

* Naming convention should include parameters used to obtain data and the dataset name istsef
* Of the variables defined, the following must be included 1: “firingRateHisto” (this file name will be modified in code to create an individual (not overlapping) histogram image too) 2: “unitWaveforms” 3: “rasterPlot” 