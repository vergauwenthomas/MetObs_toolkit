***************************
Using the GUI
***************************

A GUI (Graphical User Interface) is under construction that helps to build
a data template and explore your dataset.

The GUI can **only be launched as a local application**, or on a remote that has a graphical backend. This means that the **GUI can not be used in Google Colab notebooks!**


Why a GUI
==================================

Building a data / metadata template can sometimes be tricky. The GUI is intended to streamline this process with a visual application.
In addition to building a template, some basic functions are implemented as well.


How to lauch the GUI
======================
As explained above, the GUI can best be launched as a local python script or as a local JupyterNotebook.
To do that, make shure you have installed the Metobs-toolkit on you machine (See the introduction)


Lauch the GUI by running this code in a Python3 console or in a Jupyter notebook

.. code-block:: python

    import metobs_toolkit

    metobs_toolkit.launch_gui() #the GUI will launch
