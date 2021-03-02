from pycbio.hgdata.genePred import GenePredTbl
from pycbio.hgdata.rangeFinder import RangeFinder


def addExons(exonRf, gp):
    """"Adds exons to RangeFinder"""
    for exon in gp.exons:
        exonRf.add(gp.chrom, exon.start, exon.end, exon, gp.strand)

def exonRangeFinder(gpTbl):
    """" Have all the transcripts from GenePredTbl and create an index of all the exons """
    exonRf = RangeFinder()
    for gp in gpTbl:
        addExons(exonRf, gp)
    return exonRf




def sharedExon(exon1, exon2):
    return (exon1.gene.name2 != exon2.gene.name2) and ((exon1.start == exon2.start) or (exon1.end == exon2.end))

def findSharedGenes(exonRf, gp):
    geneNames = list()
    for exon in gp.exons:
        for overExon in exonRf.overlapping(gp.chrom, exon.start, exon.end, gp.strand):
            if sharedExon(exon, overExon):
                if exon.gene.name2 not in geneNames:
                    geneNames.append(exon.gene.name2)
    return geneNames


def findSharedExonGenes(exonRf, gpTbl):
    for gp in gpTbl:
        findSharedGenes(exonRf, gp)
    return

def main():
    file = open('tests/cases/case1AnnotV36.gp', 'r')
    gpTbl = GenePredTbl(file)
    exonRf = exonRangeFinder(gpTbl)
    findSharedExonGenes(exonRf, gpTbl)