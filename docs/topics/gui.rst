***************************
Using the GUI
***************************

A GUI (Graphical User Interface) is under construction that helps to build
a data template and explore your dataset. This GUI is made in a seperate package: `MetObs-GUI <https://github.com/vergauwenthomas/MetObs_GUI>`_



The GUI can **only be launched as a local application**, or on a remote that has a graphical backend. This means that the **GUI can not be used in Google Colab notebooks!**

.. warning::
   The GUI is currently under development and performance can not yet be guaranteed on all OS platforms.

Why a GUI
==================================

Building a data/metadata template can sometimes be tricky. The GUI is intended to streamline this process with a visual application.
In addition to building a template, some basic functions are implemented as well.


How to launch the GUI
======================
As explained above, the GUI can best be launched as a local python script or as a local JupyterNotebook.
To do that, make sure you have installed the **Metobs-toolkit** and the **Metobs-GUI** on your machine.

.. code-block:: console

    #install the metobs-toolkit
    pip3 install metobs-toolkit
    #install the metobs-gui (currently only on github)
    pip3 install git+https://github.com/vergauwenthomas/MetObs_GUI



Launch the GUI by running this code in a Python3 console or in a Jupyter notebook

.. code-block:: python

    import metobs_gui

    metobs_gui.launch_gui() #the GUI will launch
