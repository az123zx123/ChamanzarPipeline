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