#If you have a bunch of files, then just launch it on limited numbers of files instead of all at once, otherwise it will crash on you.
#!/usr/bin/env bash

## Collect all png files in the files array
files=(/Users/dani/Dropbox/prog/Apps_CIP/JSR_2017/daniTestMSM/examples/gambier_region/orig/res/vis/tif/*tif )
## How many should be done at once
batch=50
delayparam=20

## Read the array in batches of $batch
for (( i=0; $i<${#files[@]}; i+=$batch ))
do
    ## Convert this batch
    convert -delay $delayparam -loop 0 "${files[@]:$i:$batch}" animated.$i.gif
done

## Now, merge them into a single file
convert  animated.*.gif all.gif
rm -r animated*gif
