#!/bin/bash

for sample in `cat sample.list`
do
    ./pfam_scan/pfam_scan.py  -out ./4.identification/${sample}hmm_blast.hmmpy.out -outfmt csv -evalue 1e-5 -cpu 8 ./4.identification/${sample}hmm_blast.fa /mnt/d/13.pan_genome/0.software/
    grep "PF03195" ./4.identification/${sample}hmm_blast.hmmpy.out |awk -F, '{print $1}' |sort |uniq >./4.identification/${sample}hmm_blast.hmmpy.out.list
    #grep -f ./4.identification/${sample}hmm_blast.hmmpy.out.list ./4.identification/${sample}hmm_blast.fa  > ./4.identification/${sample}.list.fa
    seqkit grep -r -f ./4.identification/${sample}hmm_blast.hmmpy.out.list ./4.identification/${sample}hmm_blast.fa  -o ./4.identification/${sample}.list.fa
done
