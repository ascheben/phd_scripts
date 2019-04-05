#!/usr/bin/env bash
find . -xdev -type d -print0 |
  while IFS= read -d '' dir; do
    echo "$(find "$dir" -maxdepth 1 -print0 | grep -zc .) $dir"
  done |
  sort -rn |
  head -50
