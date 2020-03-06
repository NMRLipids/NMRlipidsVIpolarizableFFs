#!/bin/bash
op_deffile=$1
tprFile=$2
trrFile=$3
outputPath=$4

echo system | gmx trjconv -dump 0 -f $trrFile -o nowater_conf.gro -s $tprFile
resname=$(head -n 1 ${op_deffile} | awk '{print $2}')
lastlipid=$(grep $resname nowater_conf.gro  | tail -n 1 | grep -Eo '[0-9]+' | head -1)
firstlipid=$(grep $resname nowater_conf.gro | head -n 2 | tail -n 1 | grep -Eo '[0-9]+' | head -1)
rm nowater_conf.gro
numlipid=$((lastlipid-firstlipid+1))



echo "there are $numlipid of $resname lipids in this system"



while read -r name <&5
do 
for ((i=$firstlipid; i<=${lastlipid}; i++))
do

opname=$(echo ${name} | awk '{print $1}')

cname=$(echo ${name} | awk '{print $3}')
hname=$(echo ${name} | awk '{print $4}')
echo -e "a ${cname} | a ${hname} & ri $i \nq" | gmx make_ndx -f $tprFile -o ${outputPath}/foo${opname}_${i}.ndx

echo $cname"_"$hname"_&_r_"$i | gmx rotacf -P 2 -d -n ${outputPath}/foo${opname}_${i}.ndx -f $trrFile -s $tprFile -xvg none -o ${outputPath}/${i}_rotacf_$opname.xvg
rm ${outputPath}/foo${opname}_${i}.ndx
done



done 5< $op_deffile

echo $firstlipid $lastlipid
