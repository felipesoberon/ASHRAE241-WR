#!/usr/bin/env python3
"""Check reader percentiles and failure cases using small, known datasets."""
import csv
import importlib.util
import io
from pathlib import Path
import struct
import sys
import tempfile

source = Path(__file__).resolve().parents[1] / "analysis" / "read_paired.py"
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("read_paired", source)
reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reader)

root = Path(tempfile.mkdtemp(prefix="paired-reader-", dir=sys.argv[1]))
(root / "manifest.txt").write_text(
    "format=ECAIP01\nsimulations=4\nreplicates=1\ntargets=0.001,0.5\n"
    "layout=replicate,target,simulation\n", encoding="ascii")
header = b"ECAIP01\n" + struct.pack("<QQQ2d", 2, 4, 1, 0.001, 0.5)
values = struct.pack("<8d", 8, 0, 4, 4, 0, 0, 0, 2)
path = root / "ecai.bin"
path.write_bytes(header + values)
output = io.StringIO()
reader.read_percentiles(root, [0, 0.25, 0.5, 0.75, 1], output)
rows = list(csv.DictReader(io.StringIO(output.getvalue())))
assert [float(row["ecai_l_s_person"]) for row in rows[:5]] == [0, 3, 4, 5, 8]
assert [float(row["ecai_l_s_person"]) for row in rows[5:]] == [0, 0, 0, 0.5, 2]

for data in (header + values[:-8], header + values + b"extra", b"badmagic" + header[8:] + values,
             header + struct.pack("<8d", float("nan"), 0, 4, 4, 0, 0, 0, 2)):
    path.write_bytes(data)
    try:
        reader.read_percentiles(root, [0.96], io.StringIO())
    except ValueError:
        continue
    raise AssertionError("malformed dataset was accepted")
path.write_bytes(header + values)
(root / "manifest.txt").write_text(
    "format=ECAIP01\nsimulations=4\nreplicates=1\ntargets=0.001,0.1\n"
    "layout=replicate,target,simulation\n", encoding="ascii")
try:
    reader.read_percentiles(root, [0.96], io.StringIO())
except ValueError:
    pass
else:
    raise AssertionError("mismatched manifest was accepted")
print("paired reader checks passed")
