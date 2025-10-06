#!/usr/bin/bash
sed 's/,/\t/g' $1 | sed 's/\?/-1/g' 
