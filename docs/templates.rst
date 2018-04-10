Templates
---------
The input generations consist of two parts: chosing a template (see templates_)
for the kind of calculation to perform and adding some settings to that template. Notice
that the user can either pick a specific package template or provides only generic
keywords.

.. _templates:

JSON format
~~~~~~~~~~~
.. currentmodule:: qmflows.templates.templates

The Java script object notation **JSON** is a widely used data format. This format is used together with the :mod:`json` module to implement the mechanism to load/unload the templates using the function :func:`~qmflows.templates.templates.get_template`.

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


The templates in `JSON` format are translated to python dictionaries using the :func:`get_template`. This templates can be used by themselves in the calculations
or more keywords can be attached to them modifying the default values. For example the default values to carry a single point calculation using the **ORCA** quantum
package are: ::

   from qmflows.templates import singlepoint
   print(singlepoint.specific.orca)

   basis:
         basis: 	sto_sz
   method: 	
         functional: 	lda
          method: 	dft

We can easily replace the basis set just by updating the value in the |Settings| object that represent the single point calculation: ::

   from qmflows.templates import singlepoint
   singlepoint.specific.orca.basis.basis = 'DZP'

    print(singlepoint.specific.orca)

    basis: 	
      basis: 	DZP
    method: 	
      functional: 	lda
      method: 	dft
	  

.. autofunction:: get_template
