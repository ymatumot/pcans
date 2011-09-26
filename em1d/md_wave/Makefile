LIB_DIR = ..
include $(LIB_DIR)/Makefile_inc
TARGET = a.out
OBJS = const.o init.o main.o wk.o
.PHONY : all 
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o:
	$(FC) -L$(LIB_DIR) $(FFLAGS) -I $(LIB_DIR)/common -c $< -lem1d

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(FC) -L$(LIB_DIR) -o a.out $(FFLAGS) $(OBJS) -lem1d 

$(LIB_DIR)/libem1d.a : 
	cd $(LIB_DIR) ; make

main.o : init.o const.o wk.o $(LIB_DIR)/libem1d.a
init.o : const.o $(LIB_DIR)/libem1d.a

############## moment calculation ################
moment : $(LIB_DIR)/moment/mom.out
	$(LIB_DIR)/moment/mom.out ./dat/ `cd dat; \ls *_file10.dat`

$(LIB_DIR)/moment/mom.out :
	cd $(LIB_DIR)/moment ; make
###############################################

clean :
	rm -f $(TARGET) $(OBJS) *.mod dat/*.dat mom/*.dat psd/*.dat

