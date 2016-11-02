__author__ = 'Chaoren'

rise = 3.240 # Angstrom
twistAngle = 0.343359475357 # pi=180 degree

# slides between layers
slideLen = 3.834569844303501
slideAngle = 0.29052808703733513

sequence = "AAGTTTGTACGT" # from N to C terminus
# sequence = "AAA"

atPDBfilename = "../database/pna-at-aligned.pdb"
taPDBfilename = "../database/pna-ta-aligned.pdb"
gcPDBfilename = "../database/pna-gc-aligned.pdb"
cgPDBfilename = "../database/pna-cg-aligned.pdb"

atCSVfilename = "../database/pna-at.csv"
taCSVfilename = "../database/pna-ta.csv"
gcCSVfilename = "../database/pna-gc.csv"
cgCSVfilename = "../database/pna-cg.csv"

sequenceFilename = 'pna-' + sequence + '.pdb'

twoPointsLookupTable = {}
twoPointsLookupTable['A'] = (10, 21)
twoPointsLookupTable['T'] = (10, 20)
twoPointsLookupTable['G'] = (10, 22)
twoPointsLookupTable['C'] = (10, 19)

baseLenLookupTable = {}
baseLenLookupTable['A'] = 20
baseLenLookupTable['T'] = 19
baseLenLookupTable['C'] = 18
baseLenLookupTable['G'] = 21