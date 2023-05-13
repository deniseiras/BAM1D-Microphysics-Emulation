#!/bin/bash

# n1d (new 1d source), o1d (old 1d source)
n1d=( $(more Somef90files) ) 
o1d=( $(more Somef90files | sed 's!.*/!!') )  #sed is gattering the basename
o1ddir="/home/enver/Modelos/PHYSCS-2.1.0/model/source2.1_nilo_enver/"

echo ${n1d[0]}  ${o1d[0]}
echo ${n1d[3]}  ${o1d[3]}

x=0
next="y"
upx=${#o1d[@]}
while [ "$x" -lt "$upx" -a "$next" == "y" ]
do
 echo "vim ${n1d[x]} ${o1d[x]}"
 vim ${n1d[x]} ${o1ddir}${o1d[x]} -d

 read -p "Diego, Do you want to open next file? (y)"  next
x=`expr $x + 1`
done


# 2020  more Somef90files | sed 's!.*/!!'
# 2021  more Somef90files 
# 2022  more Somef90files | sed 's!.*/!!'
# 2023  cat Somef90files 

