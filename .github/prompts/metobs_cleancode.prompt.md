---
mode: edit
---

# Clean Code Prompt
Use this prompt when a code cleanup or review is called for. 


1. **Readability**: Write code that is easy to read and understand. Use meaningful variable and function names.
2. **Consistency**: Follow the project's coding style and conventions.
3. **Modularity**: Break down code into small, reusable, and testable functions or classes.
4. **Comments**: Add comments where necessary, but avoid over-commenting. Let the code speak for itself.
5. **Error Handling**: Handle errors gracefully and provide meaningful error messages.
6. **Testing**: Write tests to cover critical functionality and edge cases.
7. **Documentation**: Update documentation to reflect changes in the code.

Some specific choices to consider:
- docstrings are in the Numpy style. Do not refer to errors directly in the docstrings.
- do not write a Numpy style docstring for dunder-methods or methods decorated as property, but write a single line consice docstring instead.
- If a docstring is already present, check if it is still valid and change accordingly.
- If a docstring is already present, check for grammar and spelling errors.
- Check for typos whenever a new variable is declared, add #TYPO at the end of that line when a typo is detected in the name of the variable.
- Check if type hinting is defined for all methods and functions. Add it if it is missing or wrong.
- The `log_entry(func)` decorator (in from metobs_toolkit.backend_collection.loggingmodule) must be added to all functions and methods, add it if it is not present.
- Fix the identation of the code if needed.
- Rearrange the imports to start with python default packages, then dependency packages and lastly  local modules
- Check the spelling of inline comments, but do not remove them.



Remember, clean code is not just for machines to execute but for humans to understand and maintain.