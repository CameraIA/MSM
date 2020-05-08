# Example: mpirun -n 4 python MPI_PMRF.py -i 'data/synthetic/orig/dat/image010.dat' -s 'data/synthetic/over/dat/image_over010.dat' -o 'data/synthetic/res/image_res010.dat' -b '512,512,1' -e '10' -m '10' -l '2'

import ctypes
import sys
import os
import platform

current_path = os.path.dirname(os.path.abspath(__file__))
if platform.system() == "Darwin":
    MPI_PMRF = ctypes.CDLL(current_path + '/lib_MPI_PMRF.dylib')
else:
    MPI_PMRF = ctypes.CDLL(current_path + '/lib_MPI_PMRF.so')

LP_c_char = ctypes.POINTER(ctypes.c_char)
LP_LP_c_char = ctypes.POINTER(LP_c_char)

MPI_PMRF.open.argtypes = (ctypes.c_int, LP_LP_c_char) # argc, argv

argc = len(sys.argv)
argv = (LP_c_char * (argc + 1))()
for i, arg in enumerate(sys.argv):
    enc_arg = arg.encode('utf-8')
    argv[i] = ctypes.create_string_buffer(enc_arg)

MPI_PMRF.main(argc, argv)
