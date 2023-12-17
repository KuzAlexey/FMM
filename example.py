# EXAMPLE HOW TO USE IN PYTHON
import ctypes

script = ctypes.CDLL('../script.so')

def FMM():
    result = script.main()
    return result

FMM()
