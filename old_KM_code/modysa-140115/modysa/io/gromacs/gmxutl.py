import ctypes as C
import gmxiobase as G

def __readInt__(fd):

    byte = C.c_ubyte()

    G.LibC.read(fd, C.pointer(byte), 1)
    result = C.byte.value*16777216
    G.LibC.read(fd, C.pointer(byte), 1)
    result += byte.value*65536
    G.LibC.read(fd, C.pointer(byte), 1)
    result += byte.value*256
    G.LibC.read(fd, C.pointer(byte), 1)
    result += byte.value

    return result


def getTRNFrameSize(fileName):

    f = G.LibC.fopen(fileName, 'r')
    fd = G.LibC.fileno(f)

    entry_t = C.c_char*24
    entry = entry_t()
    G.LibC.read(fd, entry, 24)

    ir_size = __readInt__(fd)
    e_size = __readInt__(fd)
    box_size = __readInt__(fd)
    vir_size = __readInt__(fd)
    pres_size  = __readInt__(fd)
    top_size = __readInt__(fd)
    sym_size = __readInt__(fd)
    x_size = __readInt__(fd)
    v_size = __readInt__(fd)
    f_size = __readInt__(fd)
    G.LibC.fclose(f)

    return (ir_size + e_size + box_size + vir_size + pres_size +
            top_size + sym_size + x_size + v_size + f_size + 84)

