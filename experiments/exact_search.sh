#! /bin/bash

directory_name=$1
header_argument="--header"

echo "\begin{table}"
echo "\begin{center}"
echo "\begin{tabular}[!htb]{l|r|r|r}"

for i in `echo citationCiteseer.graph coPapersCiteseer.graph enron.graph loc-gowalla_edges.graph web-Google.graph`; do

    ../bin/qc --experiment=exact-search --input-file=$directory_name/$i $header_argument --latex

    if [ "$header_argument" == "--header" ]; then
        header_argument=""
    fi
done

echo "\hline"
echo "\end{tabular}"
echo "\end{center}"
echo "\caption{Run times of exact algorithms for `basename $directory_name` data set}"
echo "\end{table}"
