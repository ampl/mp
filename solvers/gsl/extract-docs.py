#!/usr/bin/env python
# This script extracts the documentation from the code.

import mmap, re

with open('amplgsl.c', 'r+b') as input, open('index.rst', 'w') as output:
  map = mmap.mmap(input.fileno(), 0)
  for i in re.finditer(r'/\*\*(.*?)\*/', map, re.DOTALL):
    s = re.sub(r'\n +\* ?', r'\n', i.group(1))
    s = re.sub(r'\$(.+?)\$', r':math:`\1`', s, re.DOTALL)
    output.write(s)
  map.close()
