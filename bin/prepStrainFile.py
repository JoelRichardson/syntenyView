
#
# prepStrainFile.py
#
# The GFF3 files from Ensemble have the MGI id embedded in the description field.
# To generate synteny blocks, need the MGI ID as the feature's ID.
#
# Also filters out features on contigs
#

import sys
import gff3
import re

# reg exp to find/capture an MGI id
mgi_re = re.compile(r'(MGI:[0-9]+)')
# 
for f in gff3.iterate(sys.stdin):
    m = mgi_re.search(f.attributes['description'])
    # keep only genes with MGI ids
    if f.type != "gene" or not m:
        continue
    # save the original ID 
    f.attributes['ensemblId']=f.ID
    # change the ID to be the MGI id
    f.ID = m.group(1)
    # exclude features on contigs.
    if len(f.seqid) > 2:
        continue
    #
    sys.stdout.write(str(f))
