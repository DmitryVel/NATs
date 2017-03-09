from __future__ import division
import re
import sys
import os, glob
from operator import itemgetter

mode = raw_input("Please choose mode (filter for filter gtf/find for find NATs): ")
entries=["exon","CDS"]
promoter_distance=1000

def main():
    if mode=="filter":
        for infile in glob.glob(os.path.join('./', '*.gff')):
            print "GTF files in the analysis directory: " + infile
        gtf_name = raw_input("Please choose GTF file: ")
        gtf_name_base=os.path.splitext(gtf_name)[0]
        gtf=open(gtf_name,"r")
        new_gtf=gtf_name_base+"_filtered.gff"
        new_gtf=filter_gtf(gtf,new_gtf)
        gtf.close()
    else:
        for infile in glob.glob(os.path.join('./', '*.gff')):
            print "GTF files in the analysis directory: " + infile
        new_gtf=raw_input("Please choose reference GTF file: ")
        main_list=find_NATs(new_gtf)
        out_name= raw_input("Please output file name: ")
        out_file=open(out_name+".txt","w")
        out_file.write("sense name"+"\t"+"chromosome"+"\t"+"overlap type"+"\t"+"antisense name"+"\t"+"overlap length"+"\t"+"overlap start"+"\t"+"overlap end"+"\n")
        for element in main_list:
            gene=element[0]
            chromo=element[1]
            entries=element[2]
            for entry in entries:
                overlap_type=entry['type']
                name=entry['name']
                length=entry['length']
                start=entry['start']
                end=entry['end']
                out_file.write(gene+"\t"+chromo+"\t"+overlap_type+"\t"+name+"\t"+str(length)+"\t"+str(start)+"\t"+str(end)+"\n")
        out_file.close()
        

def filter_gtf(gtf,new_gtf):
    old_percent=0
    N=0
    new_gtf=open(new_gtf,"w")
    print "Filtering gtf: "
    lines=gtf.readlines()
    outs=[]
    for line in lines:
        N=N+1
        percent=int(round((N/len(lines))*100))
        if old_percent!=percent:
            old_percent=percent
            print str(percent)+"%"
        split=re.split("\s+",line)
        chromo=split[0]
        entry=split[2]
        start=split[3]
        end=split[4]
        strand=split[6]
        if entry in entries:
            outs.append((line,(chromo,strand)))
    print "Sorting GTF..."
    sorted_lines=sorted(outs,key=lambda tup: (tup[1][0],tup[1][1]))
    for line in sorted_lines:
        new_gtf.write(line[0])
    new_gtf.close()
    
def find_NATs(new_gtf):
    main_list=[]
    for infile in glob.glob(os.path.join('./', '*.txt')):
        print "Text files in the analysis directory: " + infile
    inFileName = raw_input("Please enter the gene list file name: ")
    inFile=open(inFileName,"r")
    in_genes=inFile.readlines()
    gtf=open(new_gtf,"r")
    gtf_lines=gtf.readlines()
    for gene in in_genes:
        gene=gene[:-1]
        print "Looking for "+gene+" in the genome"
        out=get_position(gene,gtf_lines)
        if out!="NA":
            overlap=get_overlap(out)
            chromo=out[0][0][0]['chromosome']
            main_list.append([gene,chromo,overlap])
        else:
            print "Could not find "+gene
    return main_list

def get_position(gene,gtf_lines):
    current_chromo="none"
    found=0
    for line in gtf_lines:
        split=re.split("\s+",line)
        chromo=split[0]
        entry=split[2]
        if current_chromo!="none":
            if current_chromo!=chromo:
                if found==0:
                    print "Looking on chromosome "+chromo
                    current_chromo=chromo
                    current_chunk=[]
                else:
                    gene_chunk=get_gene(gene,current_chunk)
                    return [gene_chunk,current_chunk]
        else:
            print "Looking on chromosome "+chromo
            current_chunk=[]
            current_chromo=chromo
        start=split[3]
        end=split[4]
        strand=split[6]
        name=split[9]
        name=name[:-1]
        current_chunk.append({'chromosome':chromo,'start':start,'end':end,'strand':strand,'entry':entry,'name':name})
        if name==gene and found==0:
            print "Found "+gene+" on chromosome "+chromo
            found=1
    return "NA"

def get_gene(gene,current_chunk):
    gene_chunk=[]
    starts=[]
    ends=[]
    for element in current_chunk:
        if element['name']==gene and element['entry']=="exon":
            starts.append(int(element['start']))
            ends.append(int(element['end']))
            gene_chunk.append(element)
            strand=element['strand']
    gene_start=min(starts)
    gene_end=max(ends)
    return [gene_chunk,[gene_start,gene_end,strand]]

def if_noncoding(ref_name,chunk):
    for element in chunk:
        entry=element['entry']
        name=element['name']
        if name==ref_name and entry=="CDS":
            return False
    return True

def getOverlap(a, b):
    overlap=max(0, min(a[1], b[1]) - max(a[0], b[0]))
    start=min(a[1], b[1])
    end=max(a[0], b[0])
    return [overlap,start,end]

def get_overlap(input_list):
    outs=[]
    [gene_chunk,chunk]=input_list
    print "Looking for NATs now..."
    [gene_start,gene_end,strand]=gene_chunk[1]
    gene_chunk=gene_chunk[0]
    if strand=="+":
        promoter=[gene_start-promoter_distance,gene_start]
    else:
        promoter=[gene_end,gene_end+promoter_distance]
    for ref_entry in gene_chunk:
        ref_strand=ref_entry['strand']
        ref_start=int(ref_entry['start'])
        ref_end=int(ref_entry['end'])
        for entry in chunk:
            strand=entry['strand']
            start=int(entry['start'])
            end=int(entry['end'])
            name=entry['name']
            element_type=entry['entry']
            if strand!=ref_strand and element_type=="exon":
                if start>ref_start and start<ref_end:
                    if if_noncoding(name,chunk):
                        overlap=getOverlap((start,end), (ref_start,ref_end))
                        if overlap[0]>=50:
                            out={'type':'exonic','name':name,'length':overlap[0],'start':overlap[1],'end':overlap[2]}
                            if out not in outs:
                                outs.append(out)
                                print name+" exonic "+str(overlap[0])
                if end>ref_start and end<ref_end:
                    if if_noncoding(name,chunk):
                        overlap=getOverlap((start,end), (ref_start,ref_end))
                        if overlap[0]>=50:
                            out={'type':'exonic','name':name,'length':overlap[0],'start':overlap[1],'end':overlap[2]}
                            if out not in outs:
                                outs.append(out)
                                print name+" exonic "+str(overlap[0])
                if start>gene_start and start<gene_end:
                    if if_noncoding(name,chunk):
                        overlap=getOverlap((start,end), (gene_start,gene_end))
                        if overlap[0]>=50:
                            out={'type':'intronic','name':name,'length':overlap[0],'start':overlap[1],'end':overlap[2]}
                            if out not in outs:
                                outs.append(out)
                                print name+" intronic "+str(overlap[0])
                if end>gene_start and end<gene_end:
                    if if_noncoding(name,chunk):
                        overlap=getOverlap((start,end), (gene_start,gene_end))
                        if overlap[0]>=50:
                            out={'type':'intronic','name':name,'length':overlap[0],'start':overlap[1],'end':overlap[2]}
                            if out not in outs:
                                outs.append(out)
                                print name+" intronic "+str(overlap[0])
                if start>promoter[0] and start<promoter[1]:
                    if if_noncoding(name,chunk):
                        overlap=getOverlap((start,end), promoter)
                        if overlap[0]>=50:
                            out={'type':'promoter','name':name,'length':overlap[0],'start':overlap[1],'end':overlap[2]}
                            if out not in outs:
                                outs.append(out)
                                print name+" promoter "+str(overlap[0])
                if end>promoter[0] and end<promoter[1]:
                    if if_noncoding(name,chunk):
                        overlap=getOverlap((start,end), promoter)
                        if overlap[0]>=50:
                            out={'type':'promoter','name':name,'length':overlap[0],'start':overlap[1],'end':overlap[2]}
                            if out not in outs:
                                outs.append(out)
                                print name+" promoter "+str(overlap[0])
    return outs
                

main()
