#!/usr/bin/env python
# This script extracts the documentation from the code.

import mmap, re

with open('amplgsl.c', 'r+b') as input, open('index.rst', 'w') as output:
  map = mmap.mmap(input.fileno(), 0)
  for i in re.finditer(r'/\*\*(.*?)\*/', map, re.DOTALL):
    output.write(re.sub(r'\n +\* ?', '\n', i.group(1)))
  map.close()
