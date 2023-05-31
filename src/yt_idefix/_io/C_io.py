import re

SING_LINE_C_COMMENT_REGEXP = re.compile(r"//.*\n")
MULTILINE_C_COMMENT_REGEXP = re.compile(r"/\*[^\*/]*\*/", flags=re.DOTALL)


def strip_comments(s: str) -> str:
    s = re.sub(MULTILINE_C_COMMENT_REGEXP, "", s)
    s = re.sub(SING_LINE_C_COMMENT_REGEXP, "\n", s)
    return s
