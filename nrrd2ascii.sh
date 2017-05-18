# convert nrrd file to nrrd ascii

unu save -f nrrd -e ascii -i "$1" -o "$1"

perl formatnrrd.pl "$1" "$1"
