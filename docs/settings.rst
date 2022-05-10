Settings
--------

:class:`~qmflows.Settings` is a :class:`dict` subclass implemented in PLAMS_ and modified in *Qmflows*.
This class represents the data in a hierarchical tree-like structure. for example:

.. code:: python

    >>> from qmflows import Settings, templates

    >>> s = Settings()  # (1)

    # generic keyword
    >>> s.basis = "DZP"  #  (2)

    # "specific" allows the user to apply specific keywords for a package
    >>> s.specific.adf.basis.core = "large"  # (3)

    >>> input_settings = templates.singlepoint.overlay(s)  # (4)


The above code snippet shows how to create a :class:`~qmflows.Settings` instance object in **(1)**,
then in **(2)** the generic keyword *basis*  declares that the "DZP" should be used together with the *large* keyword
of *ADF* as shown at **(3)**.
Finally in line **(4)** the user's keywords are merged with the defaults resultin in a input like:


.. code::

    basis
        Core large
        Type DZP
    end

    integration
        accint 4.0
    end

    scf
        converge 1e-06
        iterations 100
    wnd

    xc
        lda
    end


API
~~~
.. autoclass:: qmflows.Settings
    :members:
    :special-members:
    :exclude-members: __weakref__


.. _PLAMS: https://www.scm.com/doc/plams/components/settings.html
