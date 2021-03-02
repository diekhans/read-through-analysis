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
    return (exon1.gene.name2 != exon2.gene.name2) and ((exon1.start == exon2.start) or (exon1.end == exon2.end))

def findSharedGenes(exonRf, trans):
    geneNames = list()
    for exon in trans.exons:
        for overExon in exonRf.overlapping(trans.chrom, exon.start, exon.end, trans.strand):
            if sharedExon(exon, overExon):
                if exon.gene.name2 not in geneNames:
                    geneNames.append(exon.gene.name2)
    return geneNames


def findSharedExonGenes(exonRf, transTbl):
    for trans in transTbl:
        findSharedGenes(exonRf, trans)
    return

def main():
    file = open('tests/cases/case1AnnotV36.gp', 'r')
    transTbl = GenePredTbl(file)
    exonRf = exonRangeFinder(transTbl)
    findSharedExonGenes(exonRf, transTbl)
