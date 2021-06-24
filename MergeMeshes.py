#!/usr/bin/env python 

#  Copyright (C) 2018-2021  Pedro Gomes
#  See full notice in NOTICE.md

import sys

lines = ["NDIME= "+sys.argv[1]+"\n",
         "NZONE= 2\n",
         "IZONE= 1\n"]

fid = open("mesh.su2")
fluid = fid.readlines()
fid.close()

if not fluid[-1].endswith("\n"):
  fluid[-1] += "\n"

for line in fluid:
  lines.append(line)

lines.append("IZONE= 2\n")

fid = open("solid.su2")
solid = fid.readlines()
fid.close()

for line in solid:
  lines.append(line)

fid = open("combined.su2","w")
fid.writelines(lines)
fid.close()

