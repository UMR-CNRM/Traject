Traject
=======

__*TRAcking ObJECTs*__

---

The Traject library package is a set of Python classes and functions designed to track meteorological objects (cyclones, storms, convective cells, etc) and to facilitate their manipulation. Some visualisation routines and score functions are also available.

Many tracking algorithms have already been developed for specific purposes, such as for forecasting midlatitude cyclones, for analysis of storms in satellite data, for detecting trends of tropical cyclones in reanalysis data or in climate simulations. Such algorithms apply for a specific kind of object and for specific input data. The rationale behind Traject is to propose several algorithms, that can be applied on any kind of object, and on any kind of input data. To achieve this, object definition, algorithm and input data modules are separated. Visualisation and scores are also not embedded in the track computation process. Any user who wants to add a component to Traject (object, algorithm or input format) is invited to make propositions.

Dependencies
------------

Traject is compatible with Python3.
Traject depends on the open-source library EPyGrAM (https://github.com/UMR-CNRM/EPyGrAM/) for the manipulation of meteorological model fields and formats. EPyGrAM must be installed to run Traject. Other dependancies are generic Python modules that can be installed by pip.

Installation
------------
To install Traject, put the content of the src/ directory in a directory that is declared in PYTHONPATH.

Tests
-----

Several test cases are available in the test_cases/ directory. If needed to run the tests, input data is available on demand from the code owner.

Documentation
-------------
A Traject documentation is available in pdf format on the github.

License
-------

This software is governed by the open-source [CeCILL-C](http://www.cecill.info) license under French law, cf. LICENSE.txt.
Downloading and using this code means that you have had knowledge of the CeCILL-C license and that you accept its terms.

Citation
-------

Traject is free of use for scientific research. If Traject is used for producing results for a scientific publication, the authors shall cite Traject using the information provided in CITATION.cff.

