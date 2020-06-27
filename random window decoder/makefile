NAUTY=../../yl/nauty24/
L1=-DMAXN=WORDSIZE -DWORDSIZE=64
W1=-DMAXN=WORDSIZE -DWORDSIZE=32
SUFF=W1
NAUTY_OBJS= ${NAUTY}rng.o
### replace with -O5 later on 
CXX = g++ -Wall -O2 -g -DWORDSIZE=32 -I${NAUTY} 
CC =  gcc -Wall -O2 -g -DWORDSIZE=32 -I${NAUTY}  

all: qdist dist_rnd qhpc bp code dist try gallager qbp qgall

dist_rnd: dist_rnd.C templ.H makefile
	${CXX} -DDEBUG -o dist_rnd $< -litpp -lm

qgallager: qgallager.C templ.H makefile 
	${CXX} -o qgallager $< -litpp -lm

qgall: qgall.C templ.H makefile 
	${CXX} -o qgall $< -litpp -lm

cayley: cayley.C templ.H makefile 
	${CXX} -o cayley $< -litpp -lm

gallager: gallager.C templ.H makefile 
	${CXX} -o gallager $< -litpp -lm

qbp: qbp.C templ.H makefile 
	${CXX} -o qbp $< -litpp -lm

try: try.C templ.H makefile 
	${CXX} -o try $< -litpp -lm

bp: bp.C makefile
	${CXX} -o bp $< -litpp -lm

code: code.C dist.C makefile
	${CXX} -o code $< -litpp -lm

qhpc: qhpc.C makefile
	${CXX} -o qhpc $< -litpp -lm 

mmio.o: mmio.c mmio.h makefile
	${CC} -c $< 

qdist: qdist.C mmio.o mm_gf2_read.C makefile
	${CXX} -DDEBUG -o qdist $< mmio.o mm_gf2_read.C -litpp -lm 

dist: dist.C makefile
	${CXX} -DDEBUG -o dist $< -litpp -lm 

clean: 
	rm -f *~ *.exe *.o *.tmp tmp.*

veryclean: clean
	rm -f *.out *.dat
