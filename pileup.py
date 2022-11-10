import pysam
import gzip
import sys
import os

MAX_DEPTH=1e7

bamfile = sys.argv[1]
barcodes_file = sys.argv[2]
min_barcodes = int(sys.argv[3])
sample = sys.argv[4]
outdir = sys.argv[5]

# output file locations
os.makedirs(outdir, exist_ok=True)
cov_file = gzip.open(f'{outdir}/{sample}.coverage.txt.gz','wt')
base_files = {b: gzip.open(f'{outdir}/{sample}.{b}.txt.gz','wt') for b in 'ACGT'}

# input file
with open(barcodes_file) as f:
    # f.readline() # skip first line
    depths = {l.strip() : 0 for l in f}
print(f'whitelist: {len(depths)} barcodes')

# parse barcodes
print('counting reads by barcode...')
with pysam.AlignmentFile(bamfile, "rb") as bam:
    for read in bam.fetch(region='chrM'):
        tag = read.get_tag('CB')
        if tag in depths:
            depths[tag] += 1

filt_bc = {b for b, n in depths.items() if n >= min_barcodes}
barcodes = {b:i for i, b in enumerate(filt_bc)}
n_bc = len(barcodes)

print(f'{n_bc} barcodes found with at least {min_barcodes} reads')

# Do the pileup
# len_mito = 16569

print(f'starting pileup')
bam=pysam.AlignmentFile(bamfile, "rb")
pileup = bam.pileup("chrM", max_depth=MAX_DEPTH)
for pileupcolumn in pileup:
    # if (pileupcolumn.pos < min_base) or (pileupcolumn.pos >= max_base):
    #     continue
    if (pileupcolumn.pos % 100) == 0:
        print ("Coverage at base %s = %s" %
            (pileupcolumn.pos, pileupcolumn.n))

    # create empty count dicts 
    total_depth = dict()
    fwd_depth = {b: dict() for b in 'ACTG'}
    rev_depth = {b: dict() for b in 'ACTG'}
        
    for pileupread in pileupcolumn.pileups:
        bc_tag = pileupread.alignment.get_tag("CB")
    
        # skip invalid reads
        if pileupread.is_del or pileupread.is_refskip or (bc_tag not in barcodes):
            continue
        
        # increment count dicts
        total_depth[bc_tag] = total_depth.get(bc_tag, 0) + 1
        base = pileupread.alignment.query_sequence[pileupread.query_position]
        if pileupread.alignment.is_reverse:
            rev_depth[base][bc_tag] = rev_depth[base].get(bc_tag, 0) + 1
        else:
            fwd_depth[base][bc_tag] = rev_depth[base].get(bc_tag, 0) + 1
    
    # dump lines to file
    print_pos = pileupcolumn.pos + 1 # add 1 here for 1-indexed position

    for bc, dp in total_depth.items():
        cov_file.write(f"{print_pos}\t{bc}\t{dp}\n")
    
    
    for base in 'ACTG':
        base_file = base_files[base]
        all_keys = set(rev_depth[base].keys()).union(set(fwd_depth[base].keys()))
        for bc in all_keys:
            dp_fwd = fwd_depth[base].get(bc, 0)
            dp_rev = rev_depth[base].get(bc, 0)
            base_file.write(f"{print_pos}\t{bc}\t{dp_fwd}\t{dp_rev}\n")

for f in base_files.values():
    f.close()

cov_file.close()
bam.close()
