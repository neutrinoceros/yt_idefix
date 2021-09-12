from typing import BinaryIO, Dict


def read_header(filename: str) -> str:
    with open(filename, "rb") as fh:
        return "".join(fh.readline(256).decode() for _ in range(2))


def get_field_offset_index(fh: BinaryIO) -> Dict[str, int]:  # type: ignore
    ...
