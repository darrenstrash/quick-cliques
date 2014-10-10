make clean
make

for i in  data/*/*
do 
    echo $i;
    echo "----"
    echo "Statistics:"
    echo " "
    bin/compdegen < $i > /dev/null;
    bin/printnm < $i > /dev/null;
    echo " "
    echo "Clique counts and runtimes:"
    echo " "
    bin/tomita < $i > /dev/null; 
    bin/adjlist < $i > /dev/null; 
    bin/hybrid < $i > /dev/null; 
    bin/degeneracy < $i > /dev/null; 
    echo " "
    echo " "
done
