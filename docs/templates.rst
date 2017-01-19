Templates
---------
The input generations consist of two parts: chosing a template (see templates_)
for the kind of calculation to perform and adding some settings to that template. Notice
that the user can either pick a specific package template or provides only generic
keywords.

.. _templates:

JSON format
~~~~~~~~~~~
.. currentmodule:: qmworks.templates.templates

The Java script object notation **JSON** is a widely used data format. This format is used together with the :mod:`json` module to implement the mechanism to load/unload the templates using the function :func:`~qmworks.templates.templates.get_template`.

For Example, the default parameter for a geometry optimization using ADF are given by: ::
  
	
  {
    "specific": {
        "adf": {
            "basis": {"type": "SZ"},
            "xc": {"lda": ""},
            "integration": {"accint": 4.0},
            "scf": {
        	"converge": 1e-6,
        	"iterations": 100} },
        "dftb": {
            "task": {"runtype": "SP"},
            "dftb": {"resourcesdir": "DFTB.org/3ob-3-1"} },
        "cp2k" : {
          "FORCE_EVAL": {
              "DFT": {
                  "BASIS_SET_FILE_NAME": "",
                  "MGRID": {
                      "CUTOFF": 400,
                      "NGRIDS": 4
                  },
                  "POTENTIAL_FILE_NAME": "",
                  "PRINT": {
                      "MO": {
                          "ADD_LAST"  : "NUMERIC",
                          "EACH": {
                              "QS_SCF": 0
                          },
                          "EIGENVALUES" : "",
                          "EIGENVECTORS": "",
                          "FILENAME": "./mo.data",
                          "NDIGITS": 36,
                          "OCCUPATION_NUMBERS": ""
                      }
                  },
       ................


.. autofunction:: get_template

		  
		  
.. _dictionaries:

Dictionaries
~~~~~~~~~~~~
While templates_ are as defaults, the *JSON* files stored in the dictionaries folder are use to translate the generic keywords provided by the user to specific keywords used for each package. For instance these files are used by the class :meth:`~qmworks.packages.Package`  










