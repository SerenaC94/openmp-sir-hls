#!/bin/bash

rm -r HLS_output/
rm -r panda-temp/

bambu src/accelerator.cpp \
    --top-fname=updateGridNew -v6 --compiler=I386_CLANG13 -O3 \
    --device-name=xcu55c-2Lfsvh2892-VVD --clock-period=5 -Isrc/ \
    --generate-tb=src/testbench.cpp --simulator=VERILATOR \
    --channels-number=1  --channels-type=MEM_ACC_11 --evaluation \
    --memory-allocation-policy=NO_BRAM -fopenmp -DTHREAD_NUMBER=8 |& tee log.txt