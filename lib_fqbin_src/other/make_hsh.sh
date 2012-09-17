# La entrada a fichero es el nombre del fichero sin extensión

# crea un fichero con los bloques existentes 
rm $1.hsh
zmore $1.index | awk '{if ( FNR!=1 ) print $2}' |sort -n|uniq > $1.tmp

for block in `cat $1.tmp` ; do
        minmax=`zegrep "^[^[:space:]]* $block " $1.index|awk ' \
             BEGIN { MIN="ZZZZZZZZZZZZZZZZZZZZZZZZZZZ";MAX=""} \
             {if ((MIN>$1) && ($1!="" )) MIN=$1; \
              if (MAX<$1) MAX=$1;nlines++ } \
             END {print MIN,MAX}'` 
echo $minmax $block >> $1.hsh
done
