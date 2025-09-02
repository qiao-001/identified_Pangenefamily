#!/bin/bash
for sample in `cat sample.list`
do
    hmmsearch --cut_tc --domtblout ./2.hmmsearch/${sample}.LOB.domtblout -o ./2.hmmsearch/${sample}.LOB.hmmout ./2.hmmsearch/PF0 3195.hmm ./1.database/${sample}.fa
    echo "在${sample}第一次比对"
    ####这里的参数调整blast基因的阈值
    awk '$7<1e-5 && $1 !~ /^#/ {print $0}' ./2.hmmsearch/${sample}.LOB.domtblout| awk '{print $1}'|sort -u > ./2.hmmsearch/${sample}.filter.1st
    seqkit grep -r -f ./2.hmmsearch/${sample}.filter.1st ./1.database/${sample}.fa -o ./2.hmmsearch/${sample}1st_id.fa
    clustalw -infile=./2.hmmsearch/"${sample}"1st_id.fa -output=clustal -type=PROTEIN -outfile=./2.hmmsearch/"${sample}"1st_id.aln
    echo "多序列比对!!"
    hmmbuild ./2.hmmsearch/"${sample}"_LOB.hmm ./2.hmmsearch/"${sample}"1st_id.aln
    echo "多序列比对完成"
    hmmsearch  --domtblout ./2.hmmsearch/${sample}.new_LOB.domtblout -o  ./2.hmmsearch/${sample}.new_LOB.hmmout ./2.hmmsearch/"${sample}"_LOB.hmm ./1.database/${sample}.fa

    awk '$7<1e-5 && $1 !~ /^#/ {print $0}'  ./2.hmmsearch/${sample}.new_LOB.domtblout| awk '{print $1}'|sort -u  > ./2.hmmsearch/${sample}.filter.2st
    seqkit grep -r -f ./2.hmmsearch/${sample}.filter.2st ./1.database/${sample}.fa -o ./2.hmmsearch/"${sample}"2st_id.fa
    echo "${sample}hmmsearch完成"
    ##blastɸѡ
    echo "建立balstDB"
    cd ./1.database/
    makeblastdb -in ${sample}.fa -dbtype prot -title ${sample} -out ${sample}
    cd ..
    echo "${sample}db建立完成"
    blastp -db ./1.database/${sample} -query ./3.blast/protein-matching-PF00847.fasta  -num_threads 18  -evalue 1e-5 -outfmt 6 -out ./3.blast/${sample}.Interpro_blast.txt

    awk '{print $2}' ./3.blast/${sample}.Interpro_blast.txt  |sort |uniq >./3.blast/${sample}.blast.id
    echo "获取ID"
    cat ./3.blast/${sample}.blast.id ./2.hmmsearch/${sample}.filter.2st |sort|uniq  > ./4.identification/${sample}hmm_blast.id
    echo "合并${sample}hmmsearch和blast结果"
    #wc -l ./4.identification/${sample}hmm_blast.id
    echo "${sample}合并最终序列"
    seqkit grep -r -f ./4.identification/${sample}hmm_blast.id ./1.database/${sample}.fa -o ./4.identification/${sample}hmm_blast.fa
    echo "${sample}合并序列完毕！！！，开始下一分析"
done
echo "祝贺！全部序列比对完成"