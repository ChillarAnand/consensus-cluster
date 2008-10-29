# ---- Link --------------------------- 
euclidean.so:  euclidean.o 
	gcc -shared euclidean.o -o euclidean.so  
#
#          # ---- gcc C compile ------------------
euclidean.o:  euclidean.c
	gcc -fPIC -O3 -mtune=generic -ffast-math -c euclidean.c -I/usr/include/python2.5 -I/usr/lib/python2.5/site-packages/numpy/core/include/numpy -I/usr/lib64/python2.5/site-packages/numpy/core/include/numpy

clean:
	rm euclidean.o euclidean.so

install:
	cp euclidean.so ..
