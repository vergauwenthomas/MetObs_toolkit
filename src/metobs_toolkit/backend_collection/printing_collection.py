title_char = "="
max_n_chars = 80
indent = " " * 2


def print_fmt_title(title):
    retstr = title_char * max_n_chars + "\n"
    retstr += str(title).center(max_n_chars) + "\n"
    retstr += title_char * max_n_chars + "\n"
    retstr += "\n"
    return retstr


def print_fmt_section(section):
    printsection = "--- " + str(section) + " ---"
    retstr = "\n" + printsection + "\n\n"
    return retstr


def print_fmt_line(line, identlvl=1, itemizesymbol="-"):
    """
    Print a line with indentation based on the level of indentation.
    """
    # strip leading witespace (witespace is fixed by identation)
    line = line.lstrip()
    if identlvl > 0:
        retstr = indent * identlvl + itemizesymbol + str(line)
    else:
        retstr = indent * identlvl + str(line)
    if len(retstr) >= (max_n_chars - 3):
        retstr = retstr[: max_n_chars - 3] + "..."
    retstr += "\n"
    return retstr


def print_fmt_dict(d: dict, identlvl=1, itemizesymbol="-"):
    """
    Print a dictionary with indentation based on the level of indentation.
    """
    retstr = ""
    for key, value in d.items():
        retstr += print_fmt_line(
            f"{key}: {value}", identlvl=identlvl, itemizesymbol=itemizesymbol
        )
    return retstr


def drop_title_from_str(s: str):
    """
    Drop the title from a get_info formatted string.
    """
    if s.startswith(title_char * max_n_chars):
        # drop the first 3 lines
        s = "\n".join(s.split("\n")[3:])
        return str(s)
    else:
        # string has no title
        return s
