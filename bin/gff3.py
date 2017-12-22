#
# gff3.py
#
# Utility functions for working with GFF3 files.
# This is NOT a complete implementation of GFF3! (E.g., it
# does not do validation.)
#
# Example usage: Print the ID (col 9 attribute) and length
# of all the genes in a gff file.
#
#	import gff3
#	for feature in gff3.iterate("mydata.gff"):
#	    if feature.type == "gene":
#	        length = feature.end - feature.start + 1
#	        print feature.attributes['ID'], length
#	    
# Feature objects. This module defines a class called 'Feature'. 
# A Feature object corresponds to a line in a GFF3 file and provides
# convenient access to the field values, as well as access to attribute/value
# pairs in column 9.
#
# A feature can be created several ways:
#	- With no arguments, every field is initialized to ".", except
#	col 9, which is initialized to {}.
#	- With a string argument (e.g. a line from a GFF file), the 
#	string is parsed and the fields initialized accordingly.
#	- With a Feature object, the new Feature is a copy.
#	
# Fields of a feature can be accessed:
#	- by index: f[i], where is 0 <= i <= 8
#	  Beware that indexes are 0-based
#	  while common usage in GFF documentation/discussion is 1-based.
#	  E.g., f[2] is the type column (GFF column 3).
#	- by name: f.xxx, where xxx is one of the defined column names:
#	    seqid, source, type, start, end, score, strand, phase, attributes
#	  Thus, for example, f.start is equivalent to f[3]
#
# Features also provide support for attributes in column 9.
# Attributes are maintained in a dict that maps names
# to values. A value is either a string or a list of strings.
# You can directly manipulate the attributes as a normal Python dict,
# e.g.,
#	if not f.attributes.has_key("Name"):
#	    f.attributes["Name"] = "foo"
#
# As long as an attribute name (1) does not conflict with a 
# predefined field name (above) and (2) is a valid python
# identifier, you can also access attributes as if they were
# direct attributes of the feature, e.g., f.Name is
# equivalent to f.attributes["Name"]
#
# The standard functions hasattr, getattr, and setattr work
# consistently with the above semantics.
#
#----------------------------------------------------
import sys
import types
import urllib
import re

#----------------------------------------------------
HEADER = '##gff-version 3\n'

C9SEP = ';'
QUOTECHARS_RE = re.compile(r'[\t\n\r;=%&,]')
COMMENT_CHAR = '#'
GROUPSEP = "###\n"

#----------------------------------------------------
HASH	= '#'
TAB	= '\t'
NL	= '\n'
SP	= ' '
SEMI	= ';'
EQ	= '='
COMMA	= ','

WSP_RE = re.compile(r'^\s*$')

#----------------------------------------------------
#
class ParseError(RuntimeError):
    pass

#----------------------------------------------------
#
class Feature(types.ListType):

    # These are the standard field names.
    fields = [
	"seqid",
	"source",
	"type",
	"start",
	"end",
	"score",
	"strand",
	"phase",
	"attributes",
	]
    
    # a dict that maps field names to indices
    field2index = dict(map(lambda a : (a[1],a[0]), enumerate(fields)))

    #
    def __init__(self, arg=None):
	if arg is None:
	    arg = ['.'] * 9
	    arg[-1] = []
	elif type(arg) is types.StringType:
	    arg = parse(arg)
	elif len(arg) != 9:
	    raise ValueError("Invalid initializer for GFFFeature: " \
		+ (" %d fields\n" % len(arg)) + str(arg))
	types.ListType.__init__(self,arg)
	#
	t8 = type(self[8])
	if t8 is types.StringType:
	    self[8] = parseColumn9(self[8])
	elif t8 is types.ListType:
	    self[8] = dict(self[8])
	elif t8 is types.DictType:
	    # Copy the dict.
	    # If there are list valued attrs, make sure we
	    # don't share the list obj.
	    d = {}
	    for k,v in self[8].iteritems():
		if type(v) is types.ListType:
		    d[k] = v[:]
		else:
		    d[k] = v
	    self[8] = d
	if self.start != ".":
	    self.start = int(self.start)
	if self.end != ".":
	    self.end = int(self.end)

    def __hash__(self):
	return hash(self.attributes.get('ID',None))

    def __getattr__(self, name):
	i = Feature.field2index.get(name,None)
	if i is None:
	    v = self[8].get(name,None)
	    if v is None:
		raise AttributeError(name)
	    else:
		return v
	else:
	    return self[i]

    def __setattr__(self, name, value):
	if (name=="start" or name=="end") and value != ".":
	    value = int(value)
	i = Feature.field2index.get(name,None)
	if i is None:
	    self[8][name] = value
	else:
	    self[i]=value
    
    def __str__(self):
	return format(self)

#----------------------------------------------------
# A very simple file iterator that yields a sequence
# of GFF3 Features. 
#
# Args:
#  input (file name or open file) If file name is "-", reads
#	from standard input.
#  returnGroups (boolean) If True, groups Features into lists
#	before yielding. This only makes sense if the GFF3 file
#	uses the "###" construct. (See GFF3 spec.) If False,
#	(the default), yields each Feature individually.
#
def iterate(input, returnGroups=False):
    #
    # Set up the input
    #
    closeit = False
    if type(input) is types.StringType:
	if input=="-":
	    input = sys.stdin
	else:
	    input = open(input, 'r')
	    closeit = True
    group = []
    #
    # Iterate through file.
    #
    for line in input:
	if returnGroups and line == GROUPSEP and len(group) > 0:
	    yield group
	    group = []
	elif line.startswith(COMMENT_CHAR):
	    continue
	else:
	    f = Feature(line)
	    if returnGroups:
		group.append(f)
	    else:
		yield f

    if returnGroups and len(group) > 0:
	yield group
	group = []

    #
    # Close input.
    #
    if closeit:
	input.close()

#----------------------------------------------------
def index(features):
    id2feature = {}
    for f in features:
	id = f.attributes.get("ID",None)
	if id:
	    id2feature[id] = f
    return id2feature

#----------------------------------------------------
def crossReference(features):
    id2feature = index(features)
    for f in features:
	pIds = f.attributes.get("Parent",None)
	if not pIds:
	    continue
	if type(pIds) is types.StringType:
	    pIds = [pIds]
	for pid in pIds:
	    parent = id2feature[pid]
	    f.attributes.setdefault("parent",[]).append(parent)
	    parent.attributes.setdefault("children",[]).append(f)
    return id2feature
#----------------------------------------------------
#
# Parses one line from a GFF3 file.
# Returns None if the line is a comment line. Otherwise,
# returns a list of values. If parseCol9 is True,
# (the default), the 9th column is parsed into a dict
# of name-value mappings. If False, the 9th column
# is the unparsed string.
#
def parse(line, parseCol9=True):
    if line[0:1] == HASH:
        return None
    tokens = line.split(TAB)
    if len(tokens) != 9:
        raise ParseError("Wrong number of columns (%d)\n%s" % (len(tokens),line))
    if tokens[8][-1:] == NL:
        tokens[8] = tokens[8][:-1]
    if parseCol9:
	tokens[8] = parseColumn9(tokens[8])
    return tokens

#----------------------------------------------------
#
# Parses a string of name-value attributes, as defined by GFF3. 
# Returns the corresponding dictionary. 
# 
def parseColumn9(value):
    if value == ".":
	return {}
    c9 = {}
    for t in value.split(SEMI):
	if WSP_RE.match(t):
	    continue
	tt = t.split(EQ)
	if len(tt) != 2:
	    raise ParseError("Bad column 9 format near '%s'."%t)
	n = unquote(tt[0].strip())
	v = map(unquote, tt[1].strip().split(COMMA))
	if len(v) == 1:
	    v = v[0]
	c9[n] = v
    return c9

#----------------------------------------------------
#
# Substitutes the %XX hex code for the special characters: 
# tab, newline, formfeed, ampersand, equals, semicolon,
# percent, comma.
#
def quote(v):
    return QUOTECHARS_RE.sub(lambda m:"%%%0x"%ord(m.group(0)), str(v))

#----------------------------------------------------
#
# Unquotes all hex quoted characters.
#
def unquote(v):
    return urllib.unquote(str(v)) 

#----------------------------------------------------
#
PRE = ['ID','Name','Parent','Dbxref']

#----------------------------------------------------
#
# Formats a dictionary of name/value pairs appropriately
# for column 9.
#
def formatColumn9(vals):
    if type(vals) is types.StringType:
	return quote(vals)
    parts = []
    for n in PRE:
	x = vals.get(n, None)
	if x:
	    parts.append(formatAttribute(n, x))
    for n,v in vals.iteritems():
	if n not in PRE:
	    parts.append( formatAttribute(n,v) )
    ret = C9SEP.join(parts)
    return ret

#----------------------------------------------------
#
# Formats one name/value pair appropriate for inclusion in
# column 9.
#
def formatAttribute(n, v):
    if type(v) is types.ListType:
	return "%s=%s" % (quote(n), COMMA.join(map(quote,v)))
    else:
	return "%s=%s" % (quote(n), quote(v))

#----------------------------------------------------
#
# Formats a list into a GFF3 line (newline included)
#
def format(tokens):
    lt = len(tokens)
    if lt > 9:
	tokens2 = tokens[0:9]
    elif lt < 9:
	tokens2 = tokens + ['.']*(lt-9)
	tokens2[-1] = {}
    else:
	tokens2 = tokens[:]
    tokens2[8] = formatColumn9(tokens[8])
    return TAB.join(map(str,tokens2)) + NL

#----------------------------------------------------
#
#
if __name__=="__main__":
    def printeval(expr, ns):
	v = eval(expr,ns)
	print expr, "\t=", v
	return v

    def selftest():
	f = printeval('Feature()', globals())
	f = printeval("Feature('12	MGI	gene	12345678	12347890	.	+	.	ID=MGI:222222;Name=Abc')", globals())
	ns = {'f':f}
	printeval('f[0]', ns)
	printeval('f[0:4]', ns)
	printeval('f.seqid', ns)
	printeval('f.end - f.start', ns)
	printeval('f.attributes', ns)
	printeval('f.attributes["ID"]', ns)
	printeval('f.ID', ns)
	printeval('Feature.fields', globals())
	return f

    if len(sys.argv) == 2:
	if sys.argv[1] == "-":
	    fd = sys.stdin
	else:
	    fd = open(sys.argv[1])
	for line in fd:
	    print "----------------"
	    print line,
	    if line.startswith(COMMENT_CHAR):
		continue
	    f = Feature(line)
	    print f,
	fd.close()
    else:
	f=selftest()
