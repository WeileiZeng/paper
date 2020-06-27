#ifndef UTIL_H
#define UTIL_H

#if 0
// portion from nauty.h 

#if  WORDSIZE==32
EXTDEF_CLASS
#endif

#ifdef SETWORD_LONG
#define MSK3232 0xFFFFFFFF00000000UL
#define MSK1648 0xFFFF000000000000UL
#define MSK0856 0xFF00000000000000UL
#define MSK1632 0x0000FFFF00000000UL
#define MSK0840     0xFF0000000000UL
#define MSK1616         0xFFFF0000UL 
#define MSK0824         0xFF000000UL 
#define MSK0808             0xFF00UL 
#define MSK63C  0x7FFFFFFFFFFFFFFFUL
#define MSK31C          0x7FFFFFFFUL
#define MSK15C              0x7FFFUL
#define MSK64   0xFFFFFFFFFFFFFFFFUL
#define MSK32           0xFFFFFFFFUL
#define MSK16               0xFFFFUL
#endif


#if  WORDSIZE==32
#define POPCOUNT(x) (bytecount[(x)>>24 & 0xFF] + bytecount[(x)>>16 & 0xFF] \
                        + bytecount[(x)>>8 & 0xFF] + bytecount[(x) & 0xFF])
#define FIRSTBIT(x) ((x) & MSK1616 ? ((x) & MSK0824 ? \
                     leftbit[((x)>>24) & 0xFF] : 8+leftbit[(x)>>16]) \
                    : ((x) & MSK0808 ? 16+leftbit[(x)>>8] : 24+leftbit[x]))
#define BITMASK(x)  (MSK31C >> (x))
#define ALLBITS  MSK32
#define SWCHUNK0(w) ((long)((w)>>16)&0xFFFFL)
#define SWCHUNK1(w) ((long)(w)&0xFFFFL)
#endif

#endif 

#define ERROR(fmt,...)                                                 \
  do{                                                                  \
    printf("#:[31;1m *** ERROR: " fmt " ***[0m \n",##__VA_ARGS__); \
    exit(-1);                                                          \
  }                                                                    \
  while(0)

#endif /* UTIL_H */
