#!  /bin/bash

read vertices
read edges

echo -n "p edges "
echo -n "`echo $vertices | bc` "
echo `echo $edges/2 | bc`

while read edge; do
#    echo $edge
    firstvertex=`echo $edge | sed -e "s/\([[:digit:]]\+\),\([[:digit:]]\+\)/\1/g"`
    secondvertex=`echo $edge | sed -e "s/\([[:digit:]]\+\),\([[:digit:]]\+\)/\2/g"`
    if [ "$firstvertex" -lt "$secondvertex" ]; then
        firstupbyone=`echo $firstvertex + 1 | bc`
        secondupbyone=`echo $secondvertex + 1 | bc`
        echo "e $firstupbyone $secondupbyone"
    fi
done
