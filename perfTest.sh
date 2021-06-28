#!/bin/bash

for filename in /perfTest/*.ini; do
    ./tiptopCLT.py perfTest/$filename MASTSEL/data/windpsd_mavis.fits perfTest/output$filename
    
done

