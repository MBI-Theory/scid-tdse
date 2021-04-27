goal:	scid-dx.so

dx:	scid-dx.so
	dx -memory 2048 -menubar -mdf scid-dx.mdf -program scid-dx-density.net

#       -mavx -msse2 -O3
  CCS = gcc -fPIC -rdynamic -fopenmp -g \
        -mavx -O3 \
        -fstrict-aliasing -ffast-math -fno-math-errno -funsafe-math-optimizations \
        -ffinite-math-only -fno-trapping-math -funroll-loops -Wall -Wno-missing-braces
  CCL = gcc -shared -lc -dl 

  DXINC = -I/opt/dx/include
  DXOUT = -L/opt/dx/lib_linux/ -Wl,-rpath,/opt/dx/lib_linux/ -lDXlite
  DXLIB = -L/opt/dx/lib_linux/ -Wl,-rpath,/opt/dx/lib_linux/ -lDXL

  DXBIN = /opt/dx/bin

scid-dx.so:	scid-dx.o scid-dx-load.o
	$(CCL) -o scid-dx.so scid-dx.o scid-dx-load.o $(DXLIB)

scid-dx.x:	scid-dx.o
	$(CCS) -o scid-dx.x scid-dx.o $(DXOUT)

scid-dx.o:	scid-dx.c
	$(CCS) $(CFLAGS) $(DXINC) -c scid-dx.c

scid-dx-load.o:	scid-dx-load.c
	$(CCS) $(CFLAGS) $(DXINC) -o scid-dx-load.o -c scid-dx-load.c

scid-dx-load.c:	scid-dx.mdf
	$(DXBIN)/mdf2c -m scid-dx.mdf > scid-dx-load.c

clean:
	-rm -f *.o *.so *.x
