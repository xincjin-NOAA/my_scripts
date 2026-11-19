# Standalone Hilbert-curve cross-validation grouping (from GSI 2DVar / RTMA).
# phil.f90, phil1.f90, hilbert_curve.f90 are copied unchanged from
# osw_gsi/src/gsi; kinds.f90 is a minimal stand-in for GSI's kinds/constants.
FC      = gfortran   # override with: make FC=ifx
FFLAGS ?= -O2 -g -fbacktrace
OBJ = src/kinds.o src/phil.o src/phil1.o src/hilbert_curve.o src/hilbert_cv_main.o

hilbert_cv: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

src/%.o: src/%.f90
	$(FC) $(FFLAGS) -Jsrc -c $< -o $@

src/phil.o src/phil1.o: src/kinds.o
src/hilbert_curve.o: src/phil.o src/phil1.o
src/hilbert_cv_main.o: src/kinds.o src/hilbert_curve.o

clean:
	rm -f src/*.o src/*.mod hilbert_cv
