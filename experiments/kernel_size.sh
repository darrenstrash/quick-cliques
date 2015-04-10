#! /bin/bash

directory_name=$1
header_argument="--header"

echo "\begin{table}"
echo "\begin{center}"
echo "\begin{tabular}[tb]{l|r|r|r|r|r}"

for i in `ls -1 $directory_name`; do

    ../bin/qc --experiment=kernel-size --input-file=$directory_name/$i $header_argument --latex

    if [ "$header_argument" == "--header" ]; then
        header_argument=""
    fi
done

echo "\hline"
echo "\end{tabular}"
echo "\end{center}"
echo "\label{table:kernel-sizes-`basename $directory_name`}"
echo "\caption{Kernel sizes for `basename $directory_name` data set}"
echo "\end{table}"
