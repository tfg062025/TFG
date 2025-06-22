#!/bin/bash

gcc -o SA RKF78.h RKF78.c FMD.h main_SA.c -lm;

for i in {1..1};
do ./SA;
done


