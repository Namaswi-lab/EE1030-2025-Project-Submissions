CC = gcc
CFLAGS = -O2 -lm

all: svd_full

svd_full: svd_full.c
	$(CC) svd_full.c -o svd_full $(CFLAGS)

run: svd_full
	./svd_full

clean:
	rm -f svd_full