from pycbio.hgdata.genePred import GenePredTbl
from pycbio.hgdata.rangeFinder import RangeFinder


def addExons(exonRf, trans):
    """"Adds exons to RangeFinder"""
    for exon in trans.exons:
        exonRf.add(trans.chrom, exon.start, exon.end, exon, trans.strand)

def exonRangeFinder(transTbl):
    """" Have all the transcripts from GenePredTbl and create an index of all the exons """
    exonRf = RangeFinder()
    for trans in transTbl:
        addExons(exonRf, trans)
    return exonRf


def sharedExon(exon1, exon2):
    return (exon1.start == exon2.start) or (exon1.end == exon2.end)

def findSharedGenes(exonRf, trans):
    geneNames = list()
    for exon in trans.exons:
        for overExon in exonRf.overlapping(trans.chrom, exon.start, exon.end, trans.strand):
            if sharedExon(exon, overExon):
                if overExon.gene.name2 not in geneNames:
                    geneNames.append(overExon.gene.name2)
    print(trans.name, geneNames)
    return geneNames


def findSharedExonGenes(exonRf, transTbl):
    readThroughs = []
    for trans in transTbl:
        geneNames = findSharedGenes(exonRf, trans)
        if len(geneNames) > 1:
            readThroughs.append(trans)
    return readThroughs

def main():
    file = open('tests/cases/case2AnnotV36.gp', 'r')
    transTbl = GenePredTbl(file)
    exonRf = exonRangeFinder(transTbl)
    readThroughs = findSharedExonGenes(exonRf, transTbl)
    print("readThroughts", readThroughs)
main()