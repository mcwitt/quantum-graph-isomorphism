#! /bin/bash

for f in $@; do
    grep '^[01]*$' $f > tmp && mv tmp $f
    awk -v bn=${f%.*} -v ex=${f##*.} '{print $0 > bn "-" NR "." ex}' RS='\n\n' $f
done
