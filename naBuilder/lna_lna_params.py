__author__ = 'Chaoren'

rise = 3.02 # Angstrom
twistAngle = 0.4525 # pi=180 degree

# slides between layers
slideLen = 4.655
slideAngle = 0.3497

sequence = "GG" # from N to C terminus
# sequence = "AAA"

atPDBfilename = "../database/lna-lna-at-aligned.pdb"
taPDBfilename = "../database/lna-lna-ta-aligned.pdb"
gcPDBfilename = "../database/lna-lna-gc-aligned.pdb"
cgPDBfilename = "../database/lna-lna-cg-aligned.pdb"

atCSVfilename = "../database/lna-lna-at.csv"
taCSVfilename = "../database/lna-lna-ta.csv"
gcCSVfilename = "../database/lna-lna-gc.csv"
cgCSVfilename = "../database/lna-lna-cg.csv"

sequenceFilename = 'lna-' + sequence + '.pdb'

twoPointsLookupTable = {}
twoPointsLookupTable['T'] = (15, 30) # first value is the index of C1' of A, second value is the index of C1' of T
twoPointsLookupTable['A'] = (7, 37) # first value is the index of C1' of T, second value is the index of C1' of A
twoPointsLookupTable['C'] = (15, 28) # first value is the index of C1' of G, second value is the index of C1' of C
twoPointsLookupTable['G'] = (4, 37) # first value is the index of C1' of C, second value is the index of C1' of G

baseLenLookupTable = {}             # number of atoms in A, T, C, G
baseLenLookupTable['A'] = 23
baseLenLookupTable['T'] = 22
baseLenLookupTable['C'] = 22
baseLenLookupTable['G'] = 24