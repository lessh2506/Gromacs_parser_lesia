import ctypes as C

Libgmx = C.cdll.LoadLibrary("libgmx.so")
LibC   = C.cdll.LoadLibrary("libc.so.6")

Matrix = (C.c_float * 3) * 3
Vector = (C.c_float * 3)
PVector = C.POINTER(Vector)
Boolean = C.c_int

try:
    Libgmx.GromacsVersion.restype = C.c_char_p
    cver = Libgmx.GromacsVersion().split()[-1]
except AttributeError: # no GromacsVersion function found
    cver = None

gmx_version = cver
