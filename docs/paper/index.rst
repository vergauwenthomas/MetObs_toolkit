###########################
Referencing
###########################


The Metobs-toolkit is published under the MIT license, you can use the software freely.
The Metobs-toolkit (v0.2.0) was published in JOSS: `Publication <https://joss.theoj.org/papers/10.21105/joss.05916#>`_.

Citing
----------

When citing Metobs-Toolkit, you can use:

*Vergauwen et al., (2024). MetObs - a Python toolkit for using non-traditional meteorological observations. Journal of Open Source Software, 9(95), 5916, https://doi.org/10.21105/joss.05916*

or using BiBTeX:


.. code-block:: bibtex

  @article{Vergauwen2024,
   doi = {10.21105/joss.05916},
   url = {https://doi.org/10.21105/joss.05916},
   year = {2024},
   publisher = {The Open Journal},
   volume = {9},
   number = {95},
   pages = {5916},
   author = {Thomas Vergauwen and Michiel Vieijra and Andrei Covaci and Amber Jacobs and Sara Top and Wout Dewettinck and Kobe Vandelanotte and Ian Hellebosch and Steven Caluwaerts},
   title = {MetObs - a Python toolkit for using non-traditional meteorological observations}, journal = {Journal of Open Source Software}
   }   


When referring to the MetObs-Toolkit software, please mention the used version.

.. code-block:: python

   import metobs_toolkit

   print(metobs_toolkit.__version__)


Publication code as an example
-------------------------------

You can find the notebook for creating the figures used in the JOSS publication as an example here:

.. toctree::
   :maxdepth: 1

   paper_figures.ipynb



About JOSS
-----------
The `Journal of Open Source Software <https://joss.theoj.org/>`_ is a developer-friendly, open-access journal for research software packages.
