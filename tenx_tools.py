
import pysam
import os
import sys
from collections import defaultdict
import numpy as np



bam_10x='phased_possorted_bam.bam'
chr= 'chr6'
start='63000852'
end='63000853'
start='63000846'
end='63000847'
#
start='120827191'
end='120827192'



def inversion_evidence(bam_10x,chr, start,end):
    """
    Finds evidence of inversion in an interval.

    Input:
    
    Output
    """
    sam_file = pysam.Samfile(bam_10x, 'rb')
    BX_codes=[]
    supp_aln=[]
    insert_size =[]
    read_ids=[]
    for read in sam_file.fetch(chr, int(start)-1,int(end)-1):
        aln_tags=dict(read.tags)
        #if read.is_secondary:
        #    print aln_tags['SA']
        if aln_tags.has_key('BX'):
            BX_codes.append(dict(read.tags)['BX'])
        read_ids.append(read.qname)
        insert_size.append(read.isize)
        if aln_tags.has_key('SA'):
            supp_aln.append(aln_tags['SA'])
            print read.qname + " "+ str(read.isize) + " " + str(read.is_secondary) + " "+aln_tags['SA']
        else:
            print read.qname + " "+ str(read.isize) + " " + str(read.is_secondary) 
    print len(BX_codes)
    print len(read_ids)
    print supp_aln
    print "Unique reads:" +  str(len(set(read_ids)))
    print "Unique barcodes:" +  str(len(set(BX_codes)))
    print insert_size
    isize_gt_1k= np.sum(np.array(insert_size)>1000)
    print "Num. reads" + "\t" + "Num. barcodes"+"\t"+ "isize >1000"+"\t"+ "Suppl"+"\n"
    print  str(len(set(read_ids)))+ "\t"+str(len(set(BX_codes)))+ "\t"+ str(isize_gt_1k)+ "\t"+str(len(supp_aln))
    
    #return insert_size


def get_BX_barcodes(bam_10x_fh, chro,start,end):
    """
    gets BX barcodes in a given genomic location
    
    """
    BX_codes=[]
    for read in bam_10x_fh.fetch(chro, int(start)-1, int(end)-1):
        ## alignment tags have BX barcode information
        aln_tags=dict(read.tags)
        if aln_tags.has_key('BX'):
            BX_codes.append(dict(read.tags)['BX'])
    return np.array(BX_codes)


def count_barcode_overlaps(bam_10x, chr_b1,start_b1,end_b1, chr_b2,start_b2,end_b2):
    """
    counts the number of overlapping barcodes between two genomic locations
    
    Input: 10X bam file, genomic location 1 and location 2

    Ouput: the number of barcode overlaps
    """
    bam_10x_fh = pysam.Samfile(bam_10x, 'rb')
    ### get BX barcodes from region 1
    BX_codes1 = get_BX_barcodes(bam_10x_fh, chr_b1, start_b1, end_b1)
    ### get BX barcodes from region 2
    BX_codes2 = get_BX_barcodes(bam_10x_fh, chr_b2, start_b2, end_b2)
    BX_olaps = np.intersect1d(BX_codes1,BX_codes2)
    return BX_olaps.size[0]
    
