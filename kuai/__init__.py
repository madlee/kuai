import ctypes

RealNumber = ctypes.c_float
Index = ctypes.c_uint16


from os.path import join as join_path, abspath, dirname
import sys

KUAI_BIN_FOLDER = abspath(join_path(dirname(__file__), "../bin32w"))
sys.path.append(KUAI_BIN_FOLDER)

