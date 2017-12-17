
import sys
import gff3

for line in sys.stdin:
  toks = line[:-1].split("\t")
  f = gff3.Feature()
  f.ID = toks[0]
  f.source = "MouseMine"
  f.type = "gene"
  f.seqid = toks[2]
  f.start = toks[3]
  f.end = toks[4]
  f.strand = "+" if toks[5] == "+1" else "-" if toks[5] == "-1" else "."
  sys.stdout.write(str(f))

