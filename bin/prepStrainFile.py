
#
# prepStrainFile.py
#
# The GFF3 files from Ensemble have the MGI id embedded in the description field.
# To generate synteny blocks, need the MGI ID as the feature's ID.
#
# Makes the MGI id the feature ID, and saves the 
# Filters out features on contigs
#

import sys
import gff3
import re

# keep track of MGI ids we've seen to avoid dups.
seen = set()

# reg exp to find/capture an MGI id
mgi_re = re.compile(r'(MGI:[0-9]+)')
# 
for i,f in enumerate(gff3.iterate(sys.stdin)):
    # exclude features on contigs.
    if len(f.seqid) > 2:
        continue
    # the MGI id
    m = mgi_re.search(f.attributes.get('description',''))
    if m:
        mgiid = m.group(1)
        #
        if mgiid in seen:
            sys.stderr.write("Duplicate detected (%s %s %s #%d).\n"%(f.attributes.get("Name",""), mgiid, f.ID, i))
            continue
        seen.add(mgiid)
        # no need to save the original ID since there is already a gene_id field with that value.
        # f.attributes['ensemblId']=f.ID
        # change the ID to be the MGI id
        f.ID = mgiid
        #
    sys.stdout.write(str(f))
