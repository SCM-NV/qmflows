Extracting numerical properties from output files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Quantum packages simulations generate output file in different formats. For examples the SCM_ simulation suite
(:class:`~qmflows.packages.SCM.ADF` and :class:`~qmflows.packages.SCM.DFTB` in *QMFlows*) generate binary outputs,
while other packages like :class:`~qmflows.packages.cp2k_package.CP2K`, :class:`~qmflows.packages.gamess.GAMESS`
and :class:`~qmflows.packages.orca.ORCA` generate ascii text files.

*QMFlows* abstract away all the different commmunication protocols with the different output formats, allowing the user to
extract the desire property by using the convention:

.. code:: python

    >>> result = job.property

where job is the simulation perform with a given package and property is the numerical value of interest (scalar or array).


The  *QMFlows* implementation of the aforemention mechanism search in the YAML files located at ``qmflows/data/dictionaries/``
for instructions about how to read that given property from the output file. Nevertheless, *Not all the properties for a given
pacakge are implemented*. **If the property of your interest is not available you can request it in the Qmflows issues page**.


Parsers
~~~~~~~
Internally *QMFlows* uses different mechanism to extract different properties from the output files. In the case of the :class:`~qmflows.packages.SCM.ADF` and
:class:`~qmflows.packages.SCM.DFTB` packages, *QMFlows* take advantages of the python interface to kftools_ files developed by the SCM_. In the case of *XML* output,
*QMFlows* direcltly uses the python built-in xml_ reader. For the output files in text format *Qmflows* uses a mixuture of awk_ and
*parsers*.

Parsers are a robust alternative to regular expressions, parsers are modular and reusable, while
|regex| tends to be abstruse and difficult to reuse. A parser is a function that decomposes a string (or binary) into its syntactic components using some predefined rules or grammar.
The library  pyparsing_ offers all the functionality to parse strings, some detail explanation about the library can be found at docs_.




.. _pyparsing: https://pyparsing.wikispaces.com/

.. _docs: https://pythonhosted.org/pyparsing/

.. |regex| replace:: :mod:`re`

.. _SCM: https://www.scm.com/

.. _KF: https://www.scm.com/doc/Scripting/Commandline_Tools/KF_command_line_utilities.html

.. _xml: https://docs.python.org/3.5/library/xml.etree.elementtree.html

.. _awk: https://www.gnu.org/software/gawk/manual/gawk.html

.. _properties: https://github.com/SCM-NV/qmflows/tree/master/qmflows/data/dictionaries

.. _kftools: https://www.scm.com/doc/plams/scm.html#kf-files
