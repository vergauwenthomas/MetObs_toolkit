# Contributing

All contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

You can contribute in many ways:


## Types of Contributions

### Feature Requests

If you are interested in a new Feature, you can add it as an [issue](https://github.com/vergauwenthomas/MetObs_toolkit/issues)
with the https://github.com/vergauwenthomas/MetObs_toolkit/labels/enhancement label.

In the issues describe clearly the new feature or functionality you want to see implemented in the toolkit. If the request is specific to your application/data, make sure to
add a sample of your pickled Dataset (use [.save_dataset()](https://vergauwenthomas.github.io/MetObs_toolkit/_autosummary/metobs_toolkit.dataset.Dataset.html#metobs_toolkit.dataset.Dataset.save_dataset)).

Assign yourself to the issue if you want to (help) implement the new request. All help is much appreciated!

### Report Bugs
Report bugs at the as a new [issue](https://github.com/vergauwenthomas/MetObs_toolkit/issues) with the https://github.com/vergauwenthomas/MetObs_toolkit/labels/bug label.

If you are reporting a bug, please include:

* Your operating system name and version.
* The version of the MetObs-toolkit (use `metobs_toolkit.__version__`).
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.
* If possible, provide a pickled version of your dataset (use [.save_dataset()](https://vergauwenthomas.github.io/MetObs_toolkit/_autosummary/metobs_toolkit.dataset.Dataset.html#metobs_toolkit.dataset.Dataset.save_dataset)). Limit the size of the dataset as much as possible.

### Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with https://github.com/vergauwenthomas/MetObs_toolkit/labels/bug and https://github.com/vergauwenthomas/MetObs_toolkit/labels/invalid is open to whoever wants to implement it.
If you find yourself not so familiar with Python you can start by filtering to the https://github.com/vergauwenthomas/MetObs_toolkit/labels/easy and https://github.com/vergauwenthomas/MetObs_toolkit/labels/good_first_issue labels.

The https://github.com/vergauwenthomas/MetObs_toolkit/labels/tricky label indicates that this issue might affect multiple modules of the toolkit, the data structures or is technical more challenging. Contact @vergauwenthomas to discuss a *plan-of-attack* in advance.

### Implement Features

Look through the GitHub issues for features. Anything tagged with https://github.com/vergauwenthomas/MetObs_toolkit/labels/enhancement and https://github.com/vergauwenthomas/MetObs_toolkit/labels/Feature is open to whoever wants to implement it.

### Write Documentation
The MetObs-toolkit could always use more documentation or spell checkers. Use the https://github.com/vergauwenthomas/MetObs_toolkit/labels/documentation to indicate that your issue is documentation-related.

### Submit Feedback
Any form of feedback is much appreciated. The best way to send feedback is to file [an issue](https://github.com/vergauwenthomas/MetObs_toolkit/issues). If you cannot find a suitable label, you do not have to specify one.



## Get Started
Ready to make code contributions? Here is how to set up a developer's environment for the toolkit.

### Required software

The following software (or equivalent) is required to set up a developer environment:
* [Anaconda](https://anaconda.org/)
* [Pandoc](https://pandoc.org/index.html)

Make sure you have this software installed before proceeding.

### Setup a developer environment
1. Clone the MetObs-toolkit locally:

  ```
  git clone git@github.com:vergauwenthomas/MetObs_toolkit.git
  ```
2. Create a conda environment and install the required packages.
  ```
  # Setup a developers' environment
  conda create -n metobs_dev python==3.9 poetry
  conda activate metobs_dev

  #optional: install Spyder as IDE
  #conda install spyder

  # Install dependencies in the developers' environment
  cd MetObs_toolkit
  poetry install --with documentation
  ```
3. Create a branch for local development which is a copy of the **dev** branch.:
 ```
  # checkout the dev branch
  git checkout dev
  git pull

  # Create a new local branch and switch to it.
  git branch name-of-your-bugfix-or-feature
  git checkout name-of-your-bugfix-or-feature
  ```
 Now you can make local changes.
 
4. Test your changes locally. The [build_and_test.sh](https://github.com/vergauwenthomas/MetObs_toolkit/blob/master/deploiment/build_and_test.sh) script builds the package and runs a series of tests. All tests must be successful before your contributions can be merged in the dev branch.

  ```
  source deploiment/build_and_test.sh
  ```
5. Push your code online:
   ```
   # Add your changes to your commit
   git add -A
   # Write commit text
   git commit -m "Some text describing your code changes in this commit"
   # Push your branch online
   #only the first time:
   git push --set-upstream origin name-of-your-bugfix-or-feature
   #all other times
   git push
   ```

## Pull Request Guidelines
Once your branch has been *pushed* to github, you can create a *Pull request* in github. Make sure that you have **referred the corresponding issues** to the *Pull request*.
If your code adaptations are still *work-in-progress* add the https://github.com/vergauwenthomas/MetObs_toolkit/labels/WIP label to it. For each push, github will perform a list of checks (package building, version control, functionality test, os-tests, documentation build test), in order to merge your contributions these tests must all be successful.

If your code is ready for review, you can add the https://github.com/vergauwenthomas/MetObs_toolkit/labels/Ready_for_Review label to it.

After the code review, and all review marks are resolved, your contributions will be merged to the *dev* branch.

 ## Versioning/Tagging
 From time to time the *dev* branch will be merged with the master with a new [*Release tag*](https://github.com/vergauwenthomas/MetObs_toolkit/releases). The new release will be deployed to [PyPi index](https://pypi.org/project/MetObs-toolkit/) with the adequate versioning specified.

# Support
For general support or questions, you can refer them to @vergauwenthomas, or by mail to (thomas.vergauwen@meteo.be).

# Acknowledgement
This file is inspired by the [RavenPy](https://github.com/CSHS-CWRA/RavenPy) project. Thank you for the inspiration!”. 
