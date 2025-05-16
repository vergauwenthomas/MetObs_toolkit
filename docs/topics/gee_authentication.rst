.. _Using_gee:

***************************
Using Google Earth Engine
***************************

The Google Earth Engine is used to download geospatial information, and model data
to use for your dataset. This is done to avoid downloading/reprojecting/preprocessing large
geospatial datasets and to make it possible to switch easily between different datasets.

Two methods are used to download the GEE data:

* Directly to your computer --> Only for small data transfers
* To your Google Drive --> Only when the direct download is not possible.



This page will help you set up your personal Google Earth Engine authentication.
This is needed because the GEE (Google Earth Engine) can only be used if you

* Have a Google developers account (free of charge)
* Create a cloud project on your developers account (sufficient free credits for these applications)
* Enable the GEE API on your project


Here is a step-by-step guide on how to do this.

.. note::

   This guide is to obtain a basic working setup. There are a lot of ways on how to
   set up a Google cloud project, but we only cover the minimum required steps.



Setup of a Google account
==================================

If you do not have a Google account, start by creating one.



Setup of a Google developers account
=============================================================================

A Google developers account is linked to your (regular) Google account.

#. Open a browser, and log in to Google with your account.
#. Go to this website, to create a developers account: https://developers.Google.com/

   #. Click on the three vertical dots --> hit start
   #. Fill in your name and (optional) affiliations --> hit next
   #. (optional) Select your interests --> hit next
   #. (optional) Confirm newsletter subscription --> hit next


Done, you have set up a Google developer account.


Setup a cloud project on your developer account
============================================================================

You need a cloud project to make use of the Google API's. Since the data requests
by the toolkit are (in most non-commercial applications) limited, the free credits
provided by Google are sufficient.

#. Create a cloud project: https://console.cloud.Google.com/projectcreate?pli=1

   #. Choose a project name and select No organization. --> hit create
   #. (It can take a few seconds to create your project, in the "Cloud overview" you should see your project appear.)


For more information on the non-commercial policy we refer to the `Earth Engine for Noncommercial and Research Use page <https://EarthEngine.Google.com/noncommercial/>`_.



Enable API's on your project
=============================================================================
In the last step, you need to enable the use of some API's on your project.

#. Go to your project platform page: https://console.cloud.Google.com/
#. Click on "APIs & Services"
#. Click at the top on "+ ENABLE APIS AND SERVICES"

   #. Search for the 'Google Earth Engine API', click on it --> hit ENABLE
   #. Register your GEE project: https://code.EarthEngine.Google.com/register

      #. Hit "Use with a cloud project" --> hit "Unpaid usage" and select 'Academia & Research'
      #. Select "Choose an existing Google Cloud Project" --> select your project --> hit "CONTINUE TO SUMMARY"
      #. Hit "CONFIRM AND CONTINUE"



Test your GEE access
=============================================================================

.. code-block:: python

    import metobs_toolkit

    # Use the demo files, and extract LCZ from GEE

    dataset = metobs_toolkit.Dataset()
    dataset.update_file_paths(input_data_file=metobs_toolkit.demo_datafile,
                            input_metadata_file=metobs_toolkit.demo_metadatafile,
                            template_file=metobs_toolkit.demo_template)

    dataset.import_data_from_file()

    # Extract LCZ using GEE:
    dataset.get_LCZ()

    # Selecting your cloud project:
        # 1. A link will appear, click on it
        # 2. (first time only) hit 'CHOOSE PROJECT' and select your existing cloud project
        # 3. Do NOT click the read_only scopes!
        # 4. Hit 'GENERATE TOKEN' --> select your Google account --> hit 'CONTINUE'
        # 5. Select both boxes and hit 'Continue'
        # 6. An authorization code is generated, copy it.
        # 7. In your notebook, paste the code in prompted-box and hit Enter


    # The LCZ are stored in the metadf attribute of your dataset.
    print(dataset.metadf)



.. note::

   If you click on select 'read-only' scopes in the authentication, you can only
   extract small data quantities from GEE. For larger data transfers, GEE will write
   the data to file on your Google Drive, which will raise an error when you select
   'read-only' scopes.


.. warning::

   Depending on how you use the toolkit (notebook, ipython, colab, scripts, ...),
   it can happen that your credential file, which is used for authentication when
   using GEE functionality, is invalid. If that is the case, an error that looks
   like:

   `EEException: Not signed up for Earth Engine or project is not registered. Visit https://developers.Google.com/Earth-Engine/guides/access`

   is thrown.

   To update your credential file (which is saved at `~/.config/EarthEngine/credentials`),
   you can use the ``connect_to_gee()`` function and pass additional arguments.
   Here an example on how to update the credential files:

   .. code-block:: python

       import metobs_toolkit

       metobs_toolkit.connect_to_gee(force=True, #create new credentials
                                     auth_mode='notebook', # 'notebook', 'localhost', 'gcloud' (requires gcloud installed) or 'colab' (works only in colab)
                                     )



   For more information on Authentication we refer to the
   `Authentication and Initialization Guide <https://developers.Google.com/Earth-Engine/guides/auth>`_ .
