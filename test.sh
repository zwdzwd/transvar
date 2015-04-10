#!/usr/bin/env bash

IFS=
while read line; do
  if [[ $line == "$"* ]]; then
    wzresult=$(echo ${line:2} | wzcall 2>&1 | awk '(!/^\[/) && !(/^input/)' | wzwrap);
  elif [[ $line == "@"* ]]; then
    echo "$wzresult" | tee /dev/stderr;
  else
    echo "$line";
  fi;
done < README.md.in > README.md.temp
