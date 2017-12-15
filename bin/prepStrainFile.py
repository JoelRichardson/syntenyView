
import sys
import gff3
import re

mgi_re = re.compile(r'(MGI:[0-9]+)')

for f in gff3.iterate(sys.stdin):
    if f.type != "gene":
        continue
    if not "MGI:" in f.attributes.get('description',''):
        continue
    m = mgi_re.search(f.attributes['description'])
    f.attributes['ensemblId']=f.ID
    f.ID = m.group(1)
    sys.stdout.write(str(f))
