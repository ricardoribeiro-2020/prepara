FC     = gfortran
FFLAGS = -O2
TARGET = prepara
SRC    = prepara.f90

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) -o $@ $<

clean:
	rm -f $(TARGET)
