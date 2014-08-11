#!/usr/bin/env python
# This script extracts the documentation from the code.

import mmap, re

output = None
with open('../solvers/amplgsl/amplgsl.c', 'r+b') as input:
  map = mmap.mmap(input.fileno(), 0)
  for i in re.finditer(r'/\*\*(.*?)\*/', map, re.DOTALL):
    s = re.sub(r'\n +\* ?', r'\n', i.group(1))
    s = re.sub(r'\$(.+?)\$', r':math:`\1`', s, flags=re.DOTALL)
    m = re.search(r'@file (.*)', s)
    if m:
      filename = m.group(1)
      if output:
        output.close()
      output = open('amplgsl/' + filename + '.rst', 'w')
      s = s[:m.start()] + s[m.end():]
    output.write(s)
  map.close()
