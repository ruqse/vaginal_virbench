"""Byte-reproducible gzip text output.

gzip.open() stores the current time and the file name in the gzip header, so
rewriting identical content produces different bytes (and different checksums).
DeterministicGzipWriter has the same text semantics as gzip.open(path, "wt",
newline=...) -- it is the same io.TextIOWrapper, over an in-memory buffer -- and
on close writes the bytes compressed with mtime=0 and no embedded file name, so
identical content always yields identical bytes.
"""
import gzip
import io


class DeterministicGzipWriter(io.TextIOWrapper):
    """Drop-in replacement for gzip.open(path, "wt", newline=...) when writing."""

    def __init__(self, path, newline=None, encoding="utf-8"):
        self._raw_buffer = io.BytesIO()
        self._gz_path = path
        super().__init__(self._raw_buffer, encoding=encoding, newline=newline)

    def close(self):
        if not self.closed:
            self.flush()
            data = self._raw_buffer.getvalue()
            with open(self._gz_path, "wb") as raw, \
                    gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
                gz.write(data)
        super().close()
