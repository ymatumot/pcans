LIB_DIR = ..
include $(LIB_DIR)/Makefile_inc
TARGET = psd.out
OBJS = boundary.o particle.o const.o main.o fio.o
FFLAGS2 = -mcmodel=medium -shared-intel

.PHONY : all 
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o:
	$(FC) $(FFLAGS) $(FFLAGS2) -c $< 

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(FC) -o psd.out $(FFLAGS) $(FFLAGS2) $(OBJS)

main.o : const.o boundary.o fio.o particle.o

clean :
	rm -f $(TARGET) $(OBJS) *.mod

