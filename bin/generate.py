#
# generate.py
#
# Given two sets of features, one from genome A and one from genome B,
# and a file of a/b feature pairs defining corrospondence between features in A and B,
# generates synteny blocks between genome A and genome B via interpolation.
#
# 

import argparse
import gff3
import sys

class SyntenyBlockGenerator:

    def __init__ (self):
        """
        Initializes the SyntenyBlockGenerator instance.
        """
        #
        self.A = None   # list of gff3.Features
        self.B = None   # list of gff3.Features
        self.AB = None  # list of [aid,bid] pairs
        #
        self.nBlocks = 0 # number of synteny blocks created.
        #
        # Create a special object to serve as the "missing" side of an insertion/deletion block.
        #
        self.INSERTED = {
            'index' : -1,
            'ID' : '',
            'chr': '',
            'start': 0,
            'end' : 0,
            'strand' : '.',
        }

        self.initArgParser()

    def go (self):
        """
        The generator's main program. Reads the inputs, does the computation,
        and writes the synteny blocks to the output.
        """
        self.parseArgs()
        self.readFiles()
        self.prepAB()
        self.aid2feat = self.prepGff(self.A, self.a2b)
        self.bid2feat = self.prepGff(self.B, self.b2a)
        self.join()
        if self.args.debug: self.writePairs()
        self.generateBlocks()
        self.writeBlocks()

    def initArgParser (self):
        """
        Sets up the parser for the command line args.
        """
        self.parser = argparse.ArgumentParser(description='Generate synteny blocks.')
        self.parser.add_argument(
            '-A',
            required=True,
            dest="fileA",
            metavar='AFEATURES', 
            help='GFF3 file of features from genome A.')

        self.parser.add_argument(
            '-B',
            required=True,
            dest="fileB",
            metavar='BFEATURES', 
            help='GFF3 file of features from genome B.')

        self.parser.add_argument(
            '-AB',
            dest="fileAB",
            metavar='FILE', 
            help='Tab delimited, 2-column file of A/B id pairs. These pairs define correspondence between features in AFEATURES and BFEATURES. If no -AB is provided, then features correspond if they have the same ID.')

        self.parser.add_argument(
            '-d',
            dest="debug",
            action="store_true",
            default=False,
            help='Debug mode.')

    def parseArgs (self) :
        """
        """
        self.args = self.parser.parse_args()

    def readFiles (self) :
        """
        Loads the 2 GFF3 files and the AB file (if specified).
        If no AB specified, generates AB so that features with same ID correspond.
        """
        self.A = self.readGff(self.args.fileA)
        self.B = self.readGff(self.args.fileB)
        if self.args.fileAB:
            # correspondence is based on data file provided by user
            self.AB = self.readTsv(self.args.fileAB)
        else:
            # correspondence is based on shared ID.
            allIds = set([f.ID for f in self.A] + [f.ID for f in self.B])
            self.AB = [ [i,i] for i in allIds ]

    def readGff (self, fname) :
        """
        Reads a GFF3 file. Returns list of gff3.Feature objects.
        """
        return list(gff3.iterate(fname))

    def readTsv (self, fname) :
        """
        Reads a tab delimited text file.
        Returns list of records, each a list of field values.
        """
        rows = []
        fd = open(fname, 'r')
        for line in fd:
            rows.append( line[:-1].split("\t") )
        fd.close()
        return rows

    def indexAB (self) :
        """
        Creates a mapping from aid to list of corresponding bids.
        Creates a mapping from bid to list of corresponding aids.
        """
        self.a2b = {}  # map from a -> [ b's ]
        self.b2a = {}  # map from b -> [ a's ]
        for a,b in self.AB:
            self.a2b.setdefault(a,[]).append(b)
            self.b2a.setdefault(b,[]).append(a)

    def prepAB (self) :
        """
        Filters the A/B pairs to contain only the 1:1s.
        """
        # index all relationships
        self.indexAB()
        # look for the 1-1's, build new list of a,b pairs
        ab1_1 = []
        for a in self.a2b:
            if len(self.a2b[a]) == 1:
                b = self.a2b[a][0]
                if len(self.b2a[b]) == 1:
                    ab1_1.append([a,b])
        #
        self.AB = ab1_1
        # reindex with just the 1-1's
        self.indexAB()

    def prepGff (self, feats, index) :
        """
        Filters, sorts, and otherwise modifies the list of GFF3 features
        to the refined list of (feature-like) objects.
        Returns an index from ID to feature-like object.
        """
        # a. Filter for features whose ID is in the index
        n = len(feats)
        feats[:] = filter(lambda f: f.ID.startswith("MGI:") and f.ID in index, feats)
        dn_a = n - len(feats)

        # b. Sort by chr+start position.
        def gffSorter (a, b) :
            if a.seqid == b.seqid:
                return cmp(a.start, b.start)
            else:
                return cmp(a.seqid, b.seqid)
        #
        feats.sort(gffSorter)

        # c. Filter to remove any overlaps between features.
        # IS THIS IMPORTANT??
        def overlaps(a, b):
            return a.seqid == b.seqid and a.start <= b.end and a.end >= b.start
        #
        nfs = []
        pf = None
        for f in feats:
            if pf and overlaps(pf, f):
                continue
            nfs.append(f)
            pf = f
        n = len(feats)
        feats[:] = nfs
        dn_c = n - len(feats)
        
        # d. Number the features, 1, 2, 3... and project just the bits we need
        nfs = []
        for i,f in enumerate(feats) :
            nf = {
                'index'  :   i,
                'ID'     :   f.ID,
                'chr'    :   f.seqid,
                'start'  :   f.start,
                'end'    :   f.end,
                'strand' :   f.strand
            }
            nfs.append(nf)
        feats[:] = nfs;
        #
        # e. Build an index from ID to feature, and return it.
        ix = {}
        for f in feats:
            ix[f['ID']] = f
        #
        return ix

    def renumber(self):
        """
        Renumbers the features in the current list of feature pairs to fill any gaps in the
        numbering sequence.
        """
        def _renumber(which):
            self.pairs.sort(
              lambda x,y: 1 if x[which] is self.INSERTED else -1 if y[which] is self.INSERTED else cmp(x[which]['index'], y[which]['index']))
            for i,p in enumerate(self.pairs): 
                if p[which]:
                    p[which]['index'] = i
        #
        _renumber('b')
        _renumber('a')
        # leave it sorted by a


    def join (self) :
        """
        Joins the features in A to their corresponding features in B.
        Generates a list of feature pairs. 
        Update: add pairs that contain just an A or just a B (to deal with
        insertions/deletions).
        """
        self.pairs = []
        for a in self.A:
            aid = a['ID']
            bid = self.a2b.get(aid,[None])[0]
            b = self.bid2feat.get(bid, self.INSERTED)
            #
            if b is self.INSERTED: continue
            # 
            pair = {
              'a': a,
              'b': b
              }
            self.pairs.append(pair)
        #
        '''
        for b in self.B:
            bid = b['ID']
            if not bid in self.b2a:
                # b outer join row
                pair = {
                    'a': self.INSERTED,
                    'b': b
                }
        '''

        # the join step may cause genes to drop out, and it is important that the
        # sequence is unbroken for each genome
        self.renumber()

    def startBlock(self,pair):
        """
        Starts a new synteny block from the given feature pair.
        Returns the block, which is a list of 4 values:
         - Block id (integer) Block ids are assigned starting at 0.
           They have no meaning outside a given set of results.
         - Block orientation ("+" or "-") Specifies whether the A and B regions 
           of the block have the same or opposite orientations in their respective genomes.
         - Block count (integer) Records how many a/b feature pairs combined to generate this block
         - Pair (pair of feaure-like objects) Looks like this: {a:feature,b:feature}.
           Each feature (-like object) is used to record the chromosome, start, and end of the syntenic
           region in its genome.
        """
        self.nBlocks += 1
        blockId = self.nBlocks
        blockCount = 1
        ori = +1
        if pair['a'] is self.INSERTED:
            ids = set([pair['b']['ID']])
            pcopy = { 'a':self.INSERTED, 'b': pair['b'].copy() }
        elif pair['b'] is self.INSERTED:
            ids = set([pair['a']['ID']])
            pcopy = { 'b':self.INSERTED, 'a': pair['a'].copy() }
        else:
            ori = +1 if (pair['a']['strand'] == pair['b']['strand']) else -1
            ids = set([pair['a']['ID'], pair['b']['ID']])
            pcopy = pair.copy()
        return [ blockId, ori, blockCount, pcopy, ids ]

    def extendBlock(self,currPair,currBlock):
        """
        Extends the given synteny block to include the coordinate
        ranges of the given pair.
        """
        bname,ori,blockCount,pair,ids = currBlock
        currBlock[2] = blockCount+1
        if pair['a'] is not self.INSERTED:
            pair['a']['start'] = min(pair['a']['start'], currPair['a']['start'])
            pair['a']['end']   = max(pair['a']['end'],   currPair['a']['end'])
            ids.add(currPair['a']['ID'])
        if pair['b'] is not self.INSERTED:
            pair['b']['start'] = min(pair['b']['start'], currPair['b']['start'])
            pair['b']['end']   = max(pair['b']['end'],   currPair['b']['end'])
            pair['b']['index'] = currPair['b']['index']
            ids.add(currPair['b']['ID'])

    def canMerge(self,currPair,currBlock):
        """
        Returns True iff the given pair can merge with (and extend)
        the given synteny block. 
        """
        if currBlock is None:
                    return False
        bid,ori,bcount,pair,ids = currBlock
        if currPair['a'] is self.INSERTED or currPair['b'] is self.INSERTED:
            cori = +1
        else:
            cori = 1 if (currPair['a']['strand']==currPair['b']['strand']) else -1
        return \
            currPair['a']['chr'] == pair['a']['chr'] \
            and currPair['b']['chr'] == pair['b']['chr'] \
            and ori == cori \
            and (currPair['b'] is self.INSERTED or currPair['b']['index'] == pair['b']['index']+ori)

    def generateBlocks (self) :
        """
        Scans the pairs, generating synteny blocks.
        """
        self.blocks = []
        currBlock = None
        for currPair in self.pairs:
            if self.canMerge(currPair,currBlock):
                self.extendBlock(currPair,currBlock)
            else:
                currBlock = self.startBlock(currPair)
                self.blocks.append(currBlock)
            currPair['block'] = currBlock[0]

    def writePairs (self) :
        for p in self.pairs:
            a = p['a']
            b = p['b']
            r = [ a['index'], b['index'], a['ID'], a['chr'], a['start'], a['end'], a['strand'], b['ID'], b['chr'], b['start'], b['end'], b['strand'] ]
            sys.stdout.write('# ' + '\t'.join([ str(x) for x in r ]) + '\n')

    def writeBlocks(self):
        """
        Writes the blocks to stdout.
        """
        b = [
              "blockId",
              "blockCount",
              "blockOri",
              "blockRatio",
              "aChr",
              "bChr",
              "aLength",
              "bLength",
              "aStart",
              "bStart",
              "aEnd",
              "bEnd",
              "aIndex",
              "bIndex",
              "ids",
            ]
        sys.stdout.write( '\t'.join(map(lambda x:str(x),b)) + '\n' )
        for block in self.blocks:
            blkid, ori, blkcount, fields, ids = block
            alen = fields['a']['end']-fields['a']['start']+1
            blen = fields['b']['end']-fields['b']['start']+1
            blkRatio = (1.0 * min(alen,blen)) / max(alen,blen);
            r = [
              blkid,
              blkcount,
              (ori==1 and "+" or "-"),
              "%1.2f"%blkRatio,
              fields['a']['chr'],
              fields['b']['chr'],
              fields['a']['end']-fields['a']['start']+1,
              fields['b']['end']-fields['b']['start']+1,
              fields['a']['start'],
              fields['b']['start'],
              fields['a']['end'],
              fields['b']['end'],
              fields['a']['index'],
              fields['b']['index'] - (blkcount-1 if ori == 1 else 0),
              ','.join(ids),
            ]
            sys.stdout.write( '\t'.join(map(lambda x:str(x),r)) + '\n' )

#
def main () :
    sbg = SyntenyBlockGenerator()
    sbg.go()

#
main()
