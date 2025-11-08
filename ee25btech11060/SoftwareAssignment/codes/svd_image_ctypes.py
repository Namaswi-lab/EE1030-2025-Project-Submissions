# File: svd_image_ctypes.py
import ctypes, os
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

# Paths
BASE = os.path.dirname(os.path.abspath(__file__))  # codes/
ROOT = os.path.dirname(BASE)                       # SoftwareAssignment/
FIGS = os.path.join(ROOT, "figs")
TABLES = os.path.join(ROOT, "tables")

# Load shared library (adjust name for your OS if needed)
libpath = os.path.abspath(os.path.join(BASE, "libsvd.so"))
lib = ctypes.CDLL(libpath)
lib.svd_compress_wrapper.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ctypes.POINTER(ctypes.c_double)
]
lib.svd_compress_wrapper.restype = ctypes.c_int

# Load image (must be present in figs)
img_path = os.path.join(FIGS, "einstein_original.jpg")
print("Looking for image at:", img_path)
img = Image.open(img_path).convert("L")
# Optionally resize if very large; keep as-assignment: (182,186) was used earlier
img = img.resize((182, 186))  # width, height
A = np.asarray(img, dtype=np.float64)
m, n = A.shape
print("Loaded image size (rows x cols):", m, "x", n)

# Flatten input arrays for C (row-major)
A_flat = np.ascontiguousarray(A.reshape(-1).astype(np.float64))
Aout_flat = np.zeros_like(A_flat)

# List of k values
k_values = [5, 20, 50, 100]
errors = []

for k in k_values:
    print("Processing k =", k)
    ret = lib.svd_compress_wrapper(A_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                   ctypes.c_int(m), ctypes.c_int(n), ctypes.c_int(k),
                                   Aout_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    if ret != 0:
        raise RuntimeError("C function returned error code %d" % ret)

    A_k = np.clip(Aout_flat.reshape((m, n)), 0, 255)
    # Save image (PGM)
    out_path = os.path.join(FIGS, f"einstein_reconstructed_k{k}.pgm")
    Image.fromarray(A_k.astype(np.uint8)).save(out_path)
    print("Saved:", out_path)

    # Compute Frobenius norm error
    err = np.linalg.norm(A - A_k, ord='fro')
    print(f"Frobenius error for k={k}: {err:.4f}")
    errors.append((k, err))

# Save simple table values
os.makedirs(TABLES, exist_ok=True)
vals_path = os.path.join(TABLES, "table_values.txt")
with open(vals_path, "w") as f:
    for k, err in errors:
        f.write(f"{k} {err:.6f}\n")

# Also write a latex-ready table snippet
tex_path = os.path.join(TABLES, "table.tex")
with open(tex_path, "w") as f:
    f.write("\\begin{table}[h]\n\\centering\n\\begin{tabular}{|c|c|}\n\\hline\n")
    f.write("k & Frobenius Error ($\\|A - A_k\\|_F$) \\\\\\hline\n")
    for k, err in errors:
        f.write(f"{k} & {err:.4f} \\\\\\hline\n")
    f.write("\\end{tabular}\n\\caption{Approximation errors for different $k$.}\n\\end{table}\n")

print("All done. Table and images saved.")
