#! /usr/bin/env bash

for i in `ls *_*`;
    do
        NEW=`echo $i|tr '_' '-'`
        mv $i $NEW

    done
