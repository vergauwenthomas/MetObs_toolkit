============
Contributing
============

All contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Feature Requests
~~~~~~~~~~~~~~~~

If you are interested in a new Feature, you can add it as an [issue](https://github.com/vergauwenthomas/MetObs_toolkit/issues)
with the *enhanced* label.

In the issues describe clearly new feature or functionallity you want to see implemented in the toolkit. If the request is specific for your application/data, make shure to
add a sample of your data, a template and a description.

Assign yourself to the issue if you want to (help) implement the new request.


Report Bugs
~~~~~~~~~~~

Report bugs at the as a new [issue](https://github.com/vergauwenthomas/MetObs_toolkit/issues) with the *enhanced* label.

If you are reporting a bug, please include:

* Your operating system name and version.
* The version the MetObs-toolkit.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.
* If possible, provide a pickled version of your dataset (use [.save_dataset()](https://vergauwenthomas.github.io/MetObs_toolkit/_autosummary/metobs_toolkit.dataset.Dataset.html#metobs_toolkit.dataset.Dataset.save_dataset)). Limit the size of the dataset as much as possible.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "invalid" is open to whoever wants to implement it.
If you found yourself not so familiar with python you can start by filtering to the "easy"-labels.

The "tricky"-label indicates that this issue might affect multiple modules of the toolkit, the data structures or is technical more challenging. Contact @vergauwenthomas to discuss a plan-of-attack in advance.


Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement" and "feature" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

The MetObs-toolkit could always use more documentation or spell checkers. Use the "documentation"-flag to indicate that your issue is documentation related.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/CSHS-CWRA/ravenpy/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `ravenpy` for local development.

1. Fork the `ravenpy` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/ravenpy.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv ravenpy
    $ cd ravenpy/
    $ pip install -e ".[dev]"

4. To ensure a consistent style, please install the pre-commit hooks to your repo::

    $ pre-commit install

   Special style and formatting checks will be run when you commit your changes. You
   can always run the hooks on their own with::

    $ pre-commit run -a

5. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

6. When you're done making changes, check that your changes pass flake8, black, and the
   tests, including testing other Python versions with tox::

    $ flake8 ravenpy tests
    $ black --check ravenpy tests
    $ pytest tests
    $ tox

   To get flake8, black, and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. If you are editing the docs, compile and open them with::

    $ make docs
    # or to simply generate the html
    $ cd docs/
    $ make html

8. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 3.8, 3.9, 3.10, and 3.11. Check
   https://github.com/CSHS-CWRA/RavenPy/actions/workflows/main.yml
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests::

    $ pytest tests.test_ravenpy


Versioning/Tagging
------------------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in HISTORY.rst).
Then run::

    $ bumpversion patch # possible: major / minor / patch
    $ git push
    $ git push --tags

Packaging
---------

When a new version has been minted (features have been successfully integrated test coverage and stability is adequate),
maintainers should update the pip-installable package (wheel and source release) on PyPI as well as the binary on conda-forge.

The Automated Approach
~~~~~~~~~~~~~~~~~~~~~~

The simplest way to package `ravenpy` is to "publish" a version on GitHuh. GitHub CI Actions are presently configured to build the library and publish the packages on PyPI automatically.

Tagged versions will trigger a GitHub Workflow (`tag-testpypi.yml`) that will attempt to build and publish the release on `TestPyPI <https://test.pypi.org>`_.

.. note::
    Should this step fail, changes may be needed in the package; Be sure to remove this tag on GitHub and locally, address any existing problems, and recreate the tag.

To upload a new version to `PyPI <https://pypi.org/>`_, simply create a new "Published" release version on GitHub to trigger the upload workflow (`publish-pypi.yml`). When publishing on GitHub, the maintainer can either set the release notes manually (based on the `HISTORY.rst`), or set GitHub to generate release notes automatically. The choice of method is up to the maintainer.

.. warning::
    A published version on TestPyPI/PyPI can never be overwritten. Be sure to verify that the package published at https://test.pypi.org/project/ravenpy/ matches expectations before publishing a release version on GitHub.

The Manual Approach