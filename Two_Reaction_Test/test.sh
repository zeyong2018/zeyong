#!/bin/bash

for file in $(find . -name "*.c"); do
  echo Processing $file

  cat beginning.txt $file > $file.modified

  mv $file.modified $file

done
