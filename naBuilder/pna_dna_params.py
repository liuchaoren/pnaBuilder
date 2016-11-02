__author__ = 'Chaoren'


rise = 2.93080880272 # Angstrom
twistAngle = 0.507528895251 # pi=180 degree

# slides between layers
slideLen = 4.1977172541669203
slideAngle = 0.53143356088952742

sequence = "ACCCCCGGGGGT" # from 5' to 3' on DNA strand
# sequence = "GCC"

sequenceFilename = 'pna-dna-' + sequence + '.pdb'


atPDBfilename = "../database/pna-dna-at-aligned.pdb"
taPDBfilename = "../database/pna-dna-ta-aligned.pdb"
gcPDBfilename = "../database/pna-dna-gc-aligned.pdb"
cgPDBfilename = "../database/pna-dna-cg-aligned.pdb"

atCSVfilename = "../database/pna-dna-at.csv"
taCSVfilename = "../database/pna-dna-ta.csv"
gcCSVfilename = "../database/pna-dna-gc.csv"
cgCSVfilename = "../database/pna-dna-cg.csv"

sequenceFilename = 'pna-dna-' + sequence + '.pdb'

twoPointsLookupTable = {}
twoPointsLookupTable['A'] = (10, 21)
twoPointsLookupTable['T'] = (10, 20)
twoPointsLookupTable['G'] = (10, 22)
twoPointsLookupTable['C'] = (10, 19)

baseLenLookupTable = {}
baseLenLookupTable['A'] = 21
baseLenLookupTable['T'] = 20
baseLenLookupTable['C'] = 19
baseLenLookupTable['G'] = 22
