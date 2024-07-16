.. _Template api:

============
Template
============

The template class holds the information on how to interpret raw datafiles.
This is done by mapping to standards defined in the toolkit. See the :ref:`Mapping to the toolkit`_ section.


.. note::
   In pracktice, a uses does not need to interact with the `Template()` class.
   A `Template()` instance is constructed from a templatefile (json). The user
   must create this templatefile, and the `Template()` is constructed from it.

The `Template` is stored as an atribute of a `Dataset`, and it can be reached as such.

.. code-block:: python

   import metobs_toolkit

   your_dataset = metobs_toolkit.Dataset()
   your_dataset.update_settings(
    input_data_file=" ... ",
    input_metadata_file=" ... ",
    template_file=" ... ", # path to your template file (json)
    )

    # Importing the data will construct the Template
    your_dataset.import_data_from_file()

    # The Template is stored in each Dataset
    your_dataset.template



.. currentmodule:: metobs_toolkit

Constructing a templatefile
-----------------------------

.. autosummary::
   :toctree: api/

   metobs_toolkit.build_template_prompt



Common methods
--------------------

.. autosummary::
   :toctree: api/

   template.Template.get_info
