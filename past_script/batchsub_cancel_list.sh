#! /bin/bash

for i in  $(seq 26883 27152); do
    scancel $i
done