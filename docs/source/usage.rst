Usage
=====
.. _Connecting:

Connecting files
-------------------
The tutorial is based on Google colab. The first step is connecting Google Drive to a Google Colab Notebook.
For example:

.. code-block:: python

   from google.colab import drive
   drive.mount('/content/drive', force_remount=True)

Dependency
----------
.. note::

   Currently, the dependency of the pipeline is complicated and redundant.

To install and import the related libraries, both in-house and third-party, the following steps are used:

First we need to install third-party libraries:

.. code-block::

   !pip install spiketoolkit==0.7.2
   !pip install spikewidgets==0.5.1
   !pip install spikeextractors==0.8.4
   !pip install spikesorters==0.4.3
   !pip install spikecomparison==0.3.1
   !pip install spikemetrics==0.2.0
   !pip install MEAutility==1.4.6
   !pip install spikeinterface==0.11.0
   !pip install ml_ms4alg

.. note::

   The spikeinterface version in this example is outdated.

Then we import spikeinterface (0.13).

.. code-block::

   import spikeinterface as si
   import spikeinterface.extractors as se
   import spikeinterface.toolkit as st
   import spikeinterface.sorters as ss
   import spikeinterface.comparison as sc
   import spikeinterface.widgets as sw

Next we import other related libraries.

.. code-block::

   import numpy as np
   import matplotlib.pylab as plt
   import matplotlib.image as mpimg
   import pandas as pd
   import scipy.io
   from matplotlib import cm
   import os
   import sys

Finally, we import our own libraries

.. code-block::

   #setup to import our own custom libraries
   import sys
   libdir = "/content/drive/MyDrive/EPhys Recordings/lib"
   sys.path.insert(1, libdir)

   # this is how we can import our own libraries
   from jaylib.histogram import plot_firing_rate
   from jaylib.accumulate import plot_cumulative_spikes

   # dynamically re-imports module
   import importlib
   estherlib_mods = [m for m in sys.modules.keys() if 'esther' in m]
   for m in estherlib_mods:
      importlib.reload(sys.modules[m])
   from estherlib.uniform_60_filtering import uniform_60_filt
   from estherlib.returnDataInfo import *
   from estherlib.datasetObjectModule import *
   from estherlib.cachingModule import *
   from estherlib.runRoutines import *

   from stellalib.custom_raster import custom_raster
   from madilib.waveform_viewer import waveform_subplot
