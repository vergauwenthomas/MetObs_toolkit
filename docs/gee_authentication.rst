***************************
Using Google Earth Engine
***************************

The Google earth engine is used to donwload geospatial information, and modeldata
to use for your dataset. This is done to avoid donwloading/reprojecting/preprocesing large
geospatial dataset and to make it possible to swich easely between different datasets.

There are two methods that are used to download the GEE data:

* Directly to your computer --> Only for small data transfers
* To your google drive --> Only when the direct donwload is not possible.



This script will help you how to setup your personal google earth engine authentication.
This is needed because the GEE (google earth engine) can only be used, if you

* have a google developers account (free of charge)
* Create a cloud project on your developers account (sufficient free credits for these applications)
* enable the GEE API on your project


Here is a step-by-step guide on how to do this.

Note
------
this guide is to optain a basic working setup. There are a lot of ways on how
to setup a googl cloud project, we only cover the minimum required steps.



Setup of a google account
==================================

If you do not have a google account, start by creating one.



Setup of a google developpers account
=============================================================================

A google developpers account is linked to your (regular) google account.

#. open a browser, login to google with your account.
#. Go to this website, to create a developers account: https://developers.google.com/

   #. Click on the tree vertical dots --> hit start
   #. Fill in your name and (optional) affiliations --> hit next
   #. (optinal) select your intersts --> hit next
   #. (optional) confirm newsletter subscription --> hit next


Done, you have setup a google developers account


Setup a cloud project on your developers account
============================================================================

You need a cloud project to make use of the Google API's. The API's that are used by
the toolkit have quite a lot of free credentials, so you do not need to worry about
paying for these services.

#. Create a cloud project: https://console.cloud.google.com/projectcreate?pli=1

   #. Choose a project name and select No organisation. --> hit create
   #. (It can take a few seconds to create your project, in the "Cloud overview" you should see your project appear.)



Enable API's on your project
=============================================================================
In the last step you need to enable the use of some API's on you project.

#. Go to your project platform page: https://console.cloud.google.com/
#. Click on "APIs & Services"
#. Click at the top on "+ ENABLE APIS AND SERVICES"

   #. Search for the 'Google Earth Engine API', click on it --> hit ENABLE
   #. Register your GEE project: https://code.earthengine.google.com/register

      #. hit "Use with a cloud project" --> hit "Unpaid usage" and select 'Academia & Research'
      #. Select "Choose an existing Google CLoud Project" --> select your project --> hit "CONTINUE TO SUMMARY"
      #. hit "CONFIRM AND CONTINUE"



Test your GEE acces
=============================================================================

.. code-block:: python

    import metobs_toolkit

    # Use the demo files, and extract LCZ from GEE

    dataset = metobs_toolkit.Dataset()
    dataset.update_settings(input_data_file=metobs_toolkit.demo_datafile,
                            input_metadata_file=metobs_toolkit.demo_datafiles,
                            data_template_file=metobs_toolkit.demo_template,
                            metadata_template_file=metobs_toolkit.demo_template)

    dataset.import_data_from_file()

    # Extract LCZ using GEE:
    dataset.get_lcz()

    # Selecting your cloud project:
        # 1. A link will appear, click on it
        # 2. (first time only) hit 'CHOOSE PROJECT' and select your existing cloud project
        # 3. do NOT click the read_only scopes!
        # 4. hit 'GENERATE TOKEN' --> select your google account --> hit 'CONTINUE'
        # 5. Select both boxes and hit 'Continue'
        # 6. An authorization code is generated, copy it.
        # 7. In your notebook, paste the code in propted-box and hit Enter


    # The LCZ are stored in the metadf attribute of your dataset.
    print(dataset.metadf)



Note
--------
If you click on select 'read-only' scopes in the authentication, you can only
extract small data quantities from GEE. For larger data transfer, GEE will write
the data to file on your Google drive, which will raise an error when you select
'read-only' scopes.








