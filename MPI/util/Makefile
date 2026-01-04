SHELL=/bin/sh

#ifndef MPI_INST
  MPI_INST = /usr/local/mpich-3.0.4.pgi.2019.CE
#endif

LIBNAME = libparlib.a

#ifndef CC
  CC = gcc
#endif

LIB = -L$(MPI_INST)/lib -lmpich

INC = $(MPI_INST)/include

CFLAGS = -mcmodel=medium -O2 -I $(INC)

OBJS= par_assoc_buff.o \
      par_init_flag.o \
      par_init_put.o \
      par_get_float.o \
      par_get_float8.o \
      par_get_int.o \
      par_get_noblock.o \
      par_put_int.o \
      par_put_float.o \
      par_put_float8.o \
      par_retrieve_tag.o \
      par_send_noblock.o \
      par_store_tag.o \
      par_wait.o \
      zero_tag_array.o

$(LIBNAME) : $(OBJS)
	$(CC) $(CFLAGS) -c $(?:.o=.c)
	ar -rv $@ $?
	chmod 644 $@
clean:
	rm -rf *.o $(LIBNAME)
