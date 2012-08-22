#!/usr/bin/env python
# This script extracts the documentation from the code.

import mmap, re

with open('amplgsl.c', 'r+b') as input, open('index.rst', 'w') as output:
  map = mmap.mmap(input.fileno(), 0)
  for i in re.finditer(r'/\*\*(.*?)\*/', map, re.DOTALL):
    s = re.sub(r'\n +\* ?', r'\n', i.group(1))
    s = re.sub(r'\$(.+?)\$', r':math:`\1`', s, re.DOTALL)
    m = re.search(r'@file (.*)', s)
    if m:
      filename = m.group(1)
      output.close()
      output = open(filename + '.rst', 'w')
      s = s[:m.start()] + s[m.end():]
    output.write(s)
  map.close()
