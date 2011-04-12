LIB_DIR = ..
include $(LIB_DIR)/Makefile_inc
TARGET = mom.out
OBJS = boundary.o particle.o main.o mom_calc.o fio.o
FFLAGS2 = -mcmodel=medium -shared-intel

.PHONY : all 
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o:
	$(FC) $(FFLAGS) $(FFLAGS2) -c $< 

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(FC) -o mom.out $(FFLAGS) $(FFLAGS2) $(OBJS)

main.o : boundary.o fio.o mom_calc.o particle.o

clean :
	rm -f $(TARGET) $(OBJS) *.mod

