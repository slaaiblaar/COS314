CC = g++
CFLAGS = -Wall -Wextra
# Define macros
ifdef D
    CFLAGS += -DDEBUG
endif

ifdef M
    CFLAGS += -DRUN_MANY=$(M)
endif

ifdef I
    CFLAGS += -DMAX_ITERATIONS=$(I)
endif

ifdef Q
    CFLAGS += -DNOLOG
endif

ifdef L
    CFLAGS += -DSEARCH_SPACE=LARGE_SPACE
endif

# Compile and link sa.cxx
sa:
	$(CC) -o sa $(CFLAGS) sa.cxx
	./sa
	rm sa

# Compile and link ils.cxx
ils:
	$(CC) -o ils $(CFLAGS) ils.cxx
	./ils
	rm ils

sa_tldr:
	cpp -P -DRUN_MANY=1000 -DNOLOG sa.cxx -o sa_tldr.cxx
