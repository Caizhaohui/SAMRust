# Tier-0 synthetic fixtures (CI)
#
# Regenerate:
#   python scripts/prepare_fixture.py
#
# Contents:
#   small.fa[.fai]     two contigs chr1/chr2
#   small.bam[.bai]    CIGAR + flag edge cases (M/I/D/N/S/H/=/X, dup/sec/sup/QC/unmap,
#                      placed-unmapped mate, MAPQ 0/255)
#   small.cram[.crai]  same alignments via samtools view -C -T small.fa (M11)
#   small.vcf[.gz][.tbi] SNP/indel/multi-allelic + FILTER example
#   small.bcf[.csi]      same records via bcftools view -Ob
