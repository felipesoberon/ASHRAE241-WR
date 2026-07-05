# bin_reader.py
#
# Utility for reading the C++ binary raw-data format produced by
# probability_ecai --save-all (and other C++ tools that use the
# same format).
#
# Binary layout (little-endian):
#   uint32  magic = 0x41534852 ("ASHR")
#   uint32  num_categories
#   For each category:
#     uint32  name_len
#     char[name_len]  name
#     uint32  count
#     float[count]    values  (raw probability 0-1)
#
# Returns a dict {category_name: numpy array of float32}.

import numpy as np
import struct


def load_bin(path):
    """Read a C++ .bin raw-data file and return an ordered dict
    {category: np.float32 array}."""
    with open(path, "rb") as f:
        magic = struct.unpack("<I", f.read(4))[0]
        if magic != 0x41534852:
            raise ValueError(
                f"Bad magic 0x{magic:08X} in {path} "
                f"(expected 0x41534852)"
            )
        num_cat = struct.unpack("<I", f.read(4))[0]
        data = {}
        for _ in range(num_cat):
            name_len = struct.unpack("<I", f.read(4))[0]
            name = f.read(name_len).decode("utf-8")
            count = struct.unpack("<I", f.read(4))[0]
            arr = np.fromfile(f, dtype=np.float32, count=count)
            data[name] = arr
    return data


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python bin_reader.py <file.bin>")
        sys.exit(1)
    data = load_bin(sys.argv[1])
    for cat, arr in data.items():
        print(f"{cat:<20s}  N={len(arr):>8d}  "
              f"min={arr.min():.6f}  max={arr.max():.6f}")