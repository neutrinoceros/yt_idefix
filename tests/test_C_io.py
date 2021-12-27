import pytest

from yt_idefix._io import C_io


@pytest.mark.parametrize(
    "content, expected",
    (
        ("//#define GEOMETRY CARTESIAN\n", "\n"),
        (
            "#define GEOMETRY CARTESIAN\n" "/*\n" "#define GEOMETRY SPHERICAL\n" "*/\n",
            "#define GEOMETRY CARTESIAN\n\n",
        ),
        (
            "#define ONE 1\n"
            "/*\n"
            "a multiline comment\n"
            "*/\n"
            "#define TWO 2\n"
            "/*\n"
            "a second multiline comment\n"
            "*/\n"
            "#define THREE 3\n",
            "#define ONE 1\n" "\n" "#define TWO 2\n" "\n" "#define THREE 3\n",
        ),
    ),
)
def test_strip_comments(content, expected):
    res = C_io.strip_comments(content)
    assert res == expected
