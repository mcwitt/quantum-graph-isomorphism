#! /bin/bash

for f in $@; do
    grep '^[01]*$' $f > tmp && mv tmp $f
    awk '{print $0 > FILENAME "-" NR}' RS='\n\n' $f
done
