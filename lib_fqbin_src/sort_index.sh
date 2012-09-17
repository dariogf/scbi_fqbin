#!/bin/bash

# La entrada a fichero es el nombre del fichero sin extensión

mv $1 $1.bkp
echo "UMACOMPRESSEDFORMAT 1 0 0 999999999999 999999999999" > $1
zgrep -v UMACOM $1.bkp | sort -k 1 >> $1
gzip $1
mv $1.gz $1.sorted
