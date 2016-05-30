Introduction                                                                          
=============

QMWorks
-------
Research on modern computational quantum chemistry relies on a set of computational tools
to carry out calculations. The complexity of the calculations usually requires
intercommunication between the aforementioned tools, such communication is usually done
through shell scripts that try to automatize input/output actions like: launching the
computations in a cluster, reading the resulting output and feeding the relevant numerical
result to another program. Such scripts are difficult to maintain and extend, requiring a
significant programming expertise to work with them. Being then desirable a set of automatic
and extensible tools that allows to perfom complex simulations in heterogeneous hardware plaftorms.
This library tackles the construction and efficient execution of computational chemistry workflows.
This allows computational chemists to use the emerging massively parallel compute environments
in an easy manner and focus on interpretation of scientific data rather than on tedious job
submission procedures and manual data processing.

