"""Goal is to define a decorator that copies the docstring
of one funtion to another. Useful for wrapping functions as
class methods."""

import inspect
from typing import Callable, TypeVar, Any
from typing_extensions import ParamSpec

T = TypeVar("T")
P = ParamSpec("P")


def copy_doc(
    copy_func: Callable[..., Any],
    extra_param_desc: str = "",
):
    """Copies the docstring of the given function to another and optionally
    appends an extra parameter description in numpy style.

    Parameters
    ----------
    copy_func : Callable[..., Any]
        The function whose docstring will be copied.
    extra_param_desc : str, optional
        A string representing an argument description in numpy style to be
        added to the "Parameters" section of the docstring.
    Returns
    -------
    Callable
        A decorator that applies the modified docstring to the target function.

    .. code-block:: python3

        def foo():
            '''This is a foo doc string.

            Parameters
            ----------
            x : int
                An example parameter.
            '''
            ...

        @copy_doc(foo, "y : str\n    Another example parameter.")
        def bar():
            ...
    """

    def wrapped(func: Callable[P, T]) -> Callable[P, T]:
        if copy_func.__doc__:
            trgdocstr = copy_func.__doc__

            if extra_param_desc:
                trgdocstr = add_new_arg_in_docstr(trgdocstr, extra_param_desc)
            func.__doc__ = trgdocstr
        return func

    return wrapped


def add_new_arg_in_docstr(docstr, newargstr, firstline_indent=8, description_indent=4):
    newargstr = newargstr.lstrip()  # drop leading spaces
    # Remove all leading spaces from each line in the newargstr
    newargstr = "\n".join(line.lstrip() for line in newargstr.splitlines())
    # Fix identation
    # identation at the start
    newargstr = f"{' '*firstline_indent}{newargstr}"
    # identation of the description
    newargstr = newargstr.replace(
        "\n", f"\n{' '*(firstline_indent + description_indent)}"
    )

    # Find the index of the first occurrence of 'Parameters'
    parameters_index = docstr.find("Parameters")

    if parameters_index != -1:
        # Find the index of the first '\n\n' after 'Parameters'
        insert_index = docstr.find("\n\n", parameters_index)
        if insert_index != -1:
            # Insert the new argument string at the found index
            docstr = (
                docstr[:insert_index] + "\n" + newargstr + "\n" + docstr[insert_index:]
            )
        else:
            # If no '\n\n' is found, append the new argument string at the end
            docstr += "\n" + newargstr

    else:
        raise ValueError(
            f"'Parameters' not found in the docstring string: \n\n {docstr}."
        )

    return docstr


def get_function_defaults(func: callable) -> dict:
    """
    Return the keyword arguments with default values as a dictionary for
    a function.

    Parameters
    ----------
    func : callable
        The function to inspect.

    Returns
    -------
    dict
        Dictionary of keyword arguments and their default values.
    """
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }
