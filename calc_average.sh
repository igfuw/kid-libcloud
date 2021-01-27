#!/bin/sh

#FILES= " "

DIRS=()

for SEARCHFILE in output*; do
    echo $SEARCHFILE
    [[ -d $SEARCHFILE ]] && echo "found dir: ${SEARCHFILE}"
    [[ -d $SEARCHFILE ]] && DIRS+=(${SEARCHFILE})
done

echo ${DIRS[@]}

for DIR in ${DIRS[@]}; do
  echo "working with dir ${DIR}"
  cd $DIR
  FILES=""
  #for i in {1..40}
  for i in *
  do
    echo "i = $i"
    #ncdump -v surface_precipitation_rate ${i}/1D_out.nc | grep surface_precipitation_rate
    FILE="${i}/1D_out.nc"
    if [ -f $FILE ]; then
      FILES=$FILES"${FILE} "
    else
      echo "FILE $FILE does not exist"
    fi
  done
  nces -O ${FILES} 1D_out.nc
  cp 1D_out.nc ../${DIR}_average.nc
  cd ..
done
