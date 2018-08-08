LIB_DIR = .
TARGET = $(LIB_DIR)/libem1d.a
all :  $(TARGET)

$(TARGET) : 
	cd common; make

clean :
	cd common; make clean
	cd md_ebeam; make clean
	cd md_ibeam; make clean
	cd md_shock; make clean
	cd md_wave; make clean
	cd md_whistler; make clean
	rm -f $(TARGET)

