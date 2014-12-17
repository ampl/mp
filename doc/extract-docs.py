#!/usr/bin/env python
# This script extracts the documentation from the code.

import mmap, os, re

output = None
doc_dir = os.path.dirname(__file__)
amplgsl_dir = 'amplgsl'
if not os.path.exists(amplgsl_dir):
  os.mkdir(amplgsl_dir)
with open(os.path.join(doc_dir, '../src/gsl/amplgsl.c'), 'r+b') as input:
  map = mmap.mmap(input.fileno(), 0)
  for i in re.finditer(r'/\*\*(.*?)\*/', map, re.DOTALL):
    s = re.sub(r'\n +\* ?', r'\n', i.group(1))
    s = re.sub(r'\$(.+?)\$', r':math:`\1`', s, flags=re.DOTALL)
    m = re.search(r'@file (.*)', s)
    if m:
      filename = m.group(1)
      if output:
        output.close()
      output = open(os.path.join(amplgsl_dir, filename + '.rst'), 'w')
      s = s[:m.start()] + s[m.end():]
    output.write(s)
  map.close()
