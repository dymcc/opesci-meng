#/usr/bin/zsh

for FILE in $1/*.cloog 
do
    NAME=${FILE##*/}
    BASE=${NAME%.cloog}
    cloog $FILE -compilable 1 > $1/C_$BASE.c
done


