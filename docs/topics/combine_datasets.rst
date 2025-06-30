.. _combine_datasets:

***********************************************
Combining Datasets
***********************************************

When working with observational data, it is common to encounter situations where measurements from different stations or periods are stored in separate files. The MetObs toolkit expects all observations to be in a single file for creating a `Dataset`, but you can then sum these `Dataset`s to construct one (big) Dataset.

This example demonstrates how to use the sum functionality of the `Dataset` class to merge data from multiple files into a single dataset.

.. code-block:: python

    from pathlib import Path
    import metobs_toolkit

    # Define the folder containing your separate data files
    datafolder = Path('.. Path to the folder ..')  # Specify the full path to the folder containing the files
    all_csv_files = list(datafolder.glob('*.csv'))
    print(all_csv_files)

    # Create an empty Dataset
    ds = metobs_toolkit.Dataset()

    # Loop through the files and add them to the Dataset
    for csv_file in all_csv_files:
        ds_part = metobs_toolkit.Dataset()
        ds_part.import_data_from_file(
            template_file=' .. Path to template ..', 
            input_data_file=csv_file,
            input_metadata_file=' .. Path to metadatafile ..'
        )
        ds = ds + ds_part  

    # Now ds contains all observations from all files
    print(ds)


The built-in summation of the `Dataset` class is "smart". This means that the toolkit will try to concatenate
stored data, or merge data. This allows you to combine datasets that differ 
by the present stations, model data, additional metadata, present observation types, etc. 

.. warning::
    Be aware that when combining datasets, **all outliers and gaps are reset!** It is 
    highly advised to combine datasets before applying quality control or gap filling.

    (The reset is required since the estimated frequency is not guaranteed to stay
    unchanged. Since outliers and gaps depend on it, they must be reset.)


.. warning::
    In the special case where observation data (or model data) for the same observation type,
    same station and same period is present in both elements of the sum, the values
    from the right element are used! Thus a + b is not the same as b + a if and 
    only if there is an overlap in observation type, station name and timestamps.

    