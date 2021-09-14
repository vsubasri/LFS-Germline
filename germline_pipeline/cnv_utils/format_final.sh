#!/bin/bash

files=( $(find ${cnvs}/annotated -type f -name "*hg19_multianno.txt") )

for file in ${files[@]}
do
	cut -f4 ${cnvs}/final_dels > ${cnvs}/sample
	sed  -i '1i Sample' ${cnvs}/sample
	output=${cnvs}/$(basename $file .annotated.hg19_multianno.txt).final_output.txt
	paste ${cnvs}/sample $file > $output
	rm ${cnvs}/sample
done

