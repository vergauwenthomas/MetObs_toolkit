---
applyTo: "src/metobs_toolkit/*.py"
---

Format the code, and apply the following principles:

- Look for any bugs that are introduced in the codebase.

- all functions and methods (except for dunder-methods) must have a docstring. Add a suitable concise docstring in Numpy style if it is not present.
- Ensure that the docstring is up to date with the function arguments.
- if the content of a function is changed, check if the docstring is still up to date.
- Check the docstrings for grammar and spelling errors

- The `log_entry(func)` decorator (in from metobs_toolkit.backend_collection.loggingmodule) must be added to all functions and methods, add it if it is not present. If @property is present, make sure that it is applied before the @log_entry decorator.

- fix code indentation if needed

- Rearrange the imports to start with python default packages, then dependency packages and lastly  local modules

- Check the spelling of inline comments, but do not remove them.

- Check if there are clear grammar mistakes in the variable names. If so, add a line comment #TYPO to all lines where that variable is used.

- If a new class method is added, that does not start with underscore, check if it is included in the API documentation. This is done by adding it in the correct file in the docs/reference/api_reference/ folder. (ignore the docs/reference/api folder.)

- if a new class method is added, that does not start with underscore, check if a specific test for that method is present. Tests are located in the tests/ folder, using pytest framework. If not present, add it to a suitable present test, or create a new one of that class is not already tested in one of the test files. 

- apply linting 

