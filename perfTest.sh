#!/bin/bash

for pyfile in $(find perfTest -name "*.ini"); do
    pyfile1=${pyfile%.*}
    echo python tiptopCLT.py $pyfile1 $pyfile1;
    python tiptopCLT.py $pyfile1 $pyfile1;
done

