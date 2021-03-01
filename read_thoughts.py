from pycbio.hgdata.genePred import GenePredTbl
from pycbio.hgdata.rangeFinder import RangeFinder
file=open('tests/cases/case1AnnotV36.gp', 'r')
gpTbl = GenePredTbl(file)


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

exonRf=exonRangeFinder(gpTbl)


def sharedExon(exon1, exon2):
    return (exon1.gene.name2 != exon2.gene.name2) and ((exon1.start == exon2.start) or (exon1.end == exon2.end))

def overlappingExons(exonRf, gp):
    for exon in gp.exons:
        for overExon in exonRf.overlapping(gp.chrom, exon.start, exon.end, gp.strand):
            if sharedExon(exon, overExon):
                genes=list()
                genes.append(exon.gene)
                if overExon.gene not in genes:
                    transcripts=list()
                    transcripts.append(exon.gene)
                    print(transcripts)


def findSharedExonGenes(gpTbl):
    for gp in gpTbl:
        overlappingExons(exonRf, gp)
    return exonRf

exonRf=findSharedExonGenes(gpTbl)
