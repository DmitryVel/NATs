# NATs
A Python script to identify noncoding Natural Antisense Transcripts (NATs) from AceView transcriptomic database.

Usage:
1) Download AciView GTF file from:

human: ftp://ftp.ncbi.nih.gov/repository/acedb/ncbi_37_Aug10.human.genes/AceView.ncbi_37.genes_gff.gff.gz

mouse: ftp://ftp.ncbi.nih.gov/repository/acedb/ncbi_37_Sep07.mouse.genes/AceView.mm_37.genes_gff.tar.gz

2) Untar the archive into the directory with the Python script.
3) Run the script and put "filter" when asked for mode. The script will filter and sort the GTF file.
4) Prepare the list of genes of interest to look for NATs. The gene names have to be in a text format, one name per line.
5) Save the input file in the same folder with the script and the GTF file.
6) Run the script again, this type type in "find" when asked for mode.
7) Specify the reference GTF file (filtered and sorted) and the input file.
8) After search is done type in the output file name.
9) The output will be saved in the file.
