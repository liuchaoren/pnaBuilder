__author__ = 'Chaoren'
import numpy as np

pdbcolumns = (6, 11, 16, 17, 20, 22, 26, 27, 38, 46, 54, 60, 66, 76, 78, 80)
fieldnames = ("ATOM", "atom_serial_number", "atom_name", "altLoc", "resName", "chainID", "resSeq", "iCode", "x", "y",
            "z", "occupancy", "tempFactor", "segid", "element", "charge")
dtypes = {"ATOM": np.str, "atom_serial_number": np.int64, "atom_name": np.str, "altLoc": np.str, "resName": np.str, "chainID": np.str, "resSeq": np.int64, "iCode": np.str, "x": np.float64, "y": np.float64,
            "z": np.float64, "occupancy": np.float64, "tempFactor": np.float64, "segid": np.str, "element": np.str, "charge": np.float64}

pdbformatter = zip(pdbcolumns, fieldnames)

def pdb2csv(pdbfilename, csvfilename):
    pdbfile = open(pdbfilename)
    lines = pdbfile.readlines()
    csvfile = open(csvfilename, 'w')
    fieldnames_concatenate = ",".join(fieldnames)
    csvfile.write(fieldnames_concatenate + '\n')
    for eachline in lines:
        if eachline[:4] == "ATOM" or eachline[:6] == "HETATM":
            eachlineStrip = eachline[:-1] # remove \n
            fieldList = fieldSeparator(eachlineStrip)
            joinedField = ",".join([i.strip() for i in fieldList])
            csvfile.write(joinedField + '\n')
    # csvfile.write('END')
    pdbfile.close()
    csvfile.close()


def csv2pdb(dataframe, pdbfilename):
    # given a dataframe, return a file in pdb format
    pdbfile = open(pdbfilename, 'w')
    lenofdata = len(dataframe)
    dataframelocal = dataframe.fillna('')
    for i in dataframelocal.index:
        thisline = dataframelocal.ix[i]
        endcolPre = 0
        lineReconstruct = []
        for j in pdbformatter:
            endcol = j[0]
            fieldname = j[1]
            fieldwidth = endcol - endcolPre
            fieldvalue = thisline[fieldname]
            # print fieldvalue
            if fieldvalue == "ATOM":
                fieldvalue = ("%s" % fieldvalue).ljust(fieldwidth)
            elif fieldname == "atom_name":
                fieldvalue = (" " + ("%s" % fieldvalue)).ljust(fieldwidth)
            else:
                fieldvalue = ("%s" % fieldvalue).rjust(fieldwidth)
            lineReconstruct.append(fieldvalue)
            endcolPre = endcol
        lineReconstruct = ''.join(lineReconstruct) + '\n'
        pdbfile.write(lineReconstruct)
    pdbfile.write("END")
    pdbfile.close()


def fieldSeparator(line):
    # given a line of ATOM, return a list with each field being an element
    lenofline = len(line)
    fieldList = []
    endcolPre = 0
    for eachfield in pdbformatter:
        endcol = eachfield[0]
        if endcol <= lenofline:
            fieldvalue = line[endcolPre:endcol]
            fieldList.append(fieldvalue)
        else:
            fieldList.append('');
        endcolPre = endcol
    return fieldList


if __name__ == "__main__":
    line = "ATOM      1  C8' APN A   4      21.222  26.082  10.086  1.00  8.85           C"
    print fieldSeparator(line)






