# inputs_reader.py
#
# Reader for the --save-inputs binary format produced by
# probability_ecai --save-inputs.
#
# Binary layout (little-endian):
#   uint32  magic = 0x494E5054 ("INPT")
#   uint32  num_categories
#   For each category:
#     uint32  name_len
#     char[name_len]  name
#     uint32  count           (number of simulations)
#     double[count * 21]      (21 fields per sim, sequential)
#
# Field order (21 doubles per sim):
#   0:  PBR_sample    1:  lambda_bio   2:  gamma
#   3:  n_infected    4:  phi           5:  QER_sum
#   6:  mask_factor   7:  Q             8:  P
#   9:  qer.PBR_qer  10:  qer.C_drop   11:  qer.d
#  12:  qer.E        13:  qer.Vdrop    14:  qer.GVL_ml
#  15:  qer.GVL_m3   16:  qer.VF       17:  qer.RTD
#  18:  qer.DK       19:  qer.VER     20:  qer.QER_val

import numpy as np
import struct

FIELD_NAMES = [
    "PBR_sample", "lambda_bio", "gamma",
    "n_infected", "phi", "QER_sum",
    "mask_factor", "Q", "P",
    "qer_PBR_qer", "qer_C_drop", "qer_d",
    "qer_E", "qer_Vdrop", "qer_GVL_ml",
    "qer_GVL_m3", "qer_VF", "qer_RTD",
    "qer_DK", "qer_VER", "qer_QER_val",
]

N_FIELDS = len(FIELD_NAMES)  # 21


def load_inputs(path):
    """Read a --save-inputs .bin file and return a dict
    {category: 2D numpy array of shape (count, 21)}."""
    with open(path, "rb") as f:
        magic = struct.unpack("<I", f.read(4))[0]
        if magic != 0x494E5054:
            raise ValueError(
                f"Bad magic 0x{magic:08X} in {path} "
                f"(expected 0x494E5054)"
            )
        num_cat = struct.unpack("<I", f.read(4))[0]
        data = {}
        for _ in range(num_cat):
            name_len = struct.unpack("<I", f.read(4))[0]
            name = f.read(name_len).decode("utf-8")
            count = struct.unpack("<I", f.read(4))[0]
            raw = np.fromfile(f, dtype=np.float64,
                              count=count * N_FIELDS)
            data[name] = raw.reshape(count, N_FIELDS)
    return data


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python inputs_reader.py <file.bin>")
        sys.exit(1)
    data = load_inputs(sys.argv[1])
    for cat, arr in data.items():
        print(f"\n{cat}  ({arr.shape[0]} sims, {arr.shape[1]} fields)")
        # Show first 3 simulations
        for i in range(min(3, arr.shape[0])):
            print(f"  sim {i}:")
            for j, name in enumerate(FIELD_NAMES):
                print(f"    {name:<20s} = {arr[i, j]:.6g}")