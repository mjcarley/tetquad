CC = gcc -g -O0 -Wall -mavx -DHAVE_AVX_INSTRUCTIONS \
	-Werror-implicit-function-declaration -Wstrict-prototypes \
	-Wmissing-prototypes -Wmissing-declarations



CFLAGS = `gqr-config --cflags`
LIBS = `gqr-config --libs`

TQTEST = test.o tetquad.o adapt.o quadrature.o tquad.o \
	intsincos.o predicates.o

all: tqtest

tqtest: $(TQTEST)
	$(CC) $(CFLAGS) $(TQTEST) $(LIBS) -o tqtest

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.[ch]~ *.o
