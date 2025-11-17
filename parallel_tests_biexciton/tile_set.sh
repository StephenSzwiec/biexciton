#!/bin/bash


###TO DO: Need to get qq=86 correctly updated; also, if only 1 in, will not print correct
block_end=1
sub_block=1
c=1
c1=1
c2=1
subMax=1
m=0 #array counter

echo ' enter Max Number: '

read max

echo ' enter OMP thread for given job: '

read threads

#echo $max
#echo $threads

block_segment=$((max/threads))
block_end=$((block_segment + 1))

sub_block=$(((block_end - 1) / (threads - 1)))

echo 'Start of final block suggested: '
echo $block_end

echo 'Final block end:'
echo $max

echo 'Sublock segement lengths: '
echo $sub_block

#prepare to print remaining subblock count
subMax=$((threads - 1))

echo 'Loop main, equal segment block lengths'
while [ $c -le $subMax ]

do
	#echo $c

	let c1=($c-1)*$sub_block+1
	echo $c1

	let c2=$c*$sub_block
	echo $c2


	#set up records to cat in
	printf " !\$OMP SECTION\n\n  	do ae = $c1, $c2  " > rec1$c

	cat rec1$c tile_script_mid > srec$c
	
	let m=$m+1

	array[$m]='srec'$c

	let c=$c+1
done

echo 'Make final sub block: '

let c1=$c2+1
let c2=$max

printf " !\$OMP SECTION\n\n  	do ae = $c1, $c2  " > rec1$c

cat rec1$c tile_script_mid > srec$c

let m=$m+1

array[$m]='srec'$c


#may need to reset m here if a loop


echo $c1
echo $c2

#prepare concatenation
echo "cat tile_script_top" ${array[*]} "tile_script_bottom > tile_OMP_enRead.f" > cat_test_par.sh

chmod +x cat_test_par.sh

./cat_test_par.sh



rm *rec*
