Matcha
============================

Matcha is a python library to do fast barcode matching of short DNA sequences. 
It supports multi-threaded fastq input and output, and allows for arbitrary user-specified barcodes and filtering criteria.

Installation
-------------
Requirements: Python >=3.5, pybind11

Installation can be done through pip using the github source::

   pip install git+https://github.com/GreenleafLab/matcha.git


Documentation
--------------
Example usage and reference documentation are available here:


* :doc:`/example`
* :doc:`/reference`


.. toctree::
 :maxdepth: 3
 :hidden:
 
 example
 reference
 changelog

What Matcha doesn't do
-----------------------
Matcha aims to provide simple, fast, and flexible DNA barcode matching in python.
However, it currently cannot handle:
* Probabilistic error correction (e.g. using quality scores or priors about barcode abundance)
* Variable-length barcodes
* Barcodes where the read position is not known in advance
* UMI correction where the set of valid sequences is not known in advance (potentially coming soon)