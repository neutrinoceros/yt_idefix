import sys
from copy import copy


def removeprefix(s: str, prefix: str, /) -> str:
    """
    Return a str with the given prefix string removed if present.

    If the string starts with the prefix string, return string[len(prefix):].
    Otherwise, return a copy of the original string.
    """
    if sys.version_info >= (3, 9):
        return s.removeprefix(prefix)
    else:
        if s.startswith(prefix):
            return s[len(prefix) :]
        else:
            return copy(s)


def removesuffix(s: str, suffix: str, /) -> str:
    """
    Return a str with the given suffix string removed if present.

    If the string ends with the suffix string and that suffix is not empty,
    return string[:-len(suffix)]. Otherwise, return a copy of the original
    string.
    """
    if sys.version_info >= (3, 9):
        return s.removesuffix(suffix)
    else:
        if s.endswith(suffix) and suffix != "":
            return s[: -len(suffix)]
        else:
            return copy(s)
