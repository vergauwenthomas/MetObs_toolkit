.. _combine_raw_datafiles:

********************************
Combine seperate raw data files
********************************

The toolkit expects that all the observations of your dataset are stored in a single (CSV) file. This however is not always the case, an the user must combine the observations into a single file.

Here is a mini example on how this can be done using Python with `pandas` (pandas is installed as a dependecy by the MetObs toolkit, so you do not need to install it yourself).

.. code-block:: python

    from pathlib import Path #Default in python
    import pandas as pd #Dependecy in MetObs-toolkit


Reading all files one by one
------------------------------------------------
Say that we have a folder that holds 20 files, each file holds the observations of a single station. We want to read them one by one and combine them into a single dataframe.

.. code-block:: python

    # Define the folder path
    datafolder = Path('.. Path to the folder ..') #Specify the full path to the folder containing the CSV files

    # Get a list of all the CSV files in the folder
    all_csv_files = list(datafolder.glob('*.csv')) #List of all files ending with .csv in the specified folder

    print(all_csv_files)

Before we combine all the data, we must make sure that the name of the station is present
in the data chunks before we combine them. There are two common structures:

* The name of the station is present in the datafile as a seperate column. If this is the case, then we do not need to take extra actions, the combining of the data will include the names as well.

* The name of the station is not present in the datafile, but is used in the filename. In this case, we need to extract the name from the filename, add a new column in the dataframe (not manually!), and set the stationname in that column. Often the name of the station is only a part of the filename, thus we must select only the name part of the filename.


.. code-block:: python

    all_dataframes = [] #List to store all dataframes (will be filled in the loop below)

    # Loop through each CSV file and read it into a DataFrame
    for csv_file in all_csv_files:
        # Read the CSV file into a DataFrame
        df = pd.read_csv(csv_file,
                        sep=';', #CHECK THIS ! The separator used in the CSV file can differ.
                        decimal=',', #CHECK THIS ! The decimal character used in the CSV file can differ.
                        encoding='UTF-8')

        #Extract the name
        #Option 1: stationname is a column in the CSV file --> skip the following line, go to the .append(df) line!

        #Option 2: stationname is part of the file-name of the CSV file.

        filename = csv_file.stem #Extract the name of the file without the extension
        stationname = filename[-5:] #stationname are the last 5 characters of the filename

        #Create a new column 'stationname' in the DataFrame and fill it with the station name
        df['stationname'] = stationname #Add the station name to the DataFrame

        # Append the DataFrame to the list
        all_dataframes.append(df) #Add the dataframe to the list

Now we have a list (`all_dataframes`) that holds all the data as seperate dataframes. Now we combine them into a single dataframe, and write it to a csv file.



.. note::

    Do not add the combined csv file in the same folder as the seperate datafiles! If you would rerun this example, the combined-data-file will be read and merged with the other seperate data files. You will end up with a dataframe with lot's of duplicated rows.

.. code-block:: python

    fulldf = pd.concat(all_dataframes, ignore_index=True) #Concatenate all dataframes in the list into one dataframe

    # Save the concatenated DataFrame to a new CSV file
    output_file = " ... target file location ... .csv" #Specify the full path to the output CSV file
    fulldf.to_csv(output_file, sep=';', decimal=',', encoding='UTF-8', index=False) #Save the concatenated dataframe to a CSV file
    print(f"All CSV files have been concatenated and saved to {output_file}")
