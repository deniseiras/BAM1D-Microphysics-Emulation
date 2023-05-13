#!/bin/bash
num=0

for ARQ in `find . -name *.f90`
do
tmp=`wc -l | awk &acute; {print $1}&acute; `
num=`expr $num + $tmp`
done

echo Total de Linhas: $num
