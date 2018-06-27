#!/bin/bash
F90="gfortran"
SRC="npt.f90 ran1.f gasdev.f"

case $F90 in
(gfortran) OPT="-O4 -fdefault-real-8" ;;
(g95)      OPT="-O4 -r8" ;;
esac

$F90 $OPT $SRC -o npt.x
