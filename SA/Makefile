# ---- Link --------------------------- 
sa.so:  sa.o 
	gcc -shared sa.o -o sa.so  
#
#          # ---- gcc C compile ------------------
sa.o:  sa.c
	gcc -fPIC -O3 -mtune=generic -ffast-math -c sa.c -I/usr/include/python2.6 -I/usr/lib/python2.6/site-packages/numpy/core/include/numpy -I/usr/lib64/python2.6/site-packages/numpy/core/include/numpy


clean:
	rm sa.o sa.so

install:
	cp sa.so ..
