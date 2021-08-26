if [ $# -ne 4 ]
then
echo "run as haplotype assembler name samples"
fi
for s in ${@:4};
do
    samtools faidx -f ../$2_${s}/$1.scaffolds.fasta
    awk -v a=$s -v b=$3 'NR<31 {c+=$2} END {print b","a",auto_size,"c}' ../$2_${s}/$1.scaffolds.fasta.fai
    continue
    awk -v a=$s -v b=$3 '$1=="SZ" {print b","a",size,"$2}' $1_${s}_$2.contigs.auN.txt
    continue
    for i in {1..29};
    do
       awk '/Satellite/ {print $3}' ../$2_${s}/$1_split_chrm/${i}.chrm.fa.tbl; done | awk -v a=$s -v b=$3 '{c+=$1}END{print b","a",satellite,"c}'
    for i in {1..29};
    do
       awk '/bases masked/ {print $3}' ../$2_${s}/$1_split_chrm/${i}.chrm.fa.tbl; done | awk -v a=$s -v b=$3 '{c+=$1}END{print b","a",masked,"c}'
    awk -v a=$s -v b=$3 '/NG50/ {print b","a",NG50,"$2}' $1_${s}_$2.NGA50.txt
    awk '/phased/ {print $2}' $1_${s}_$2.merqury.full.stats | sed 's/,//g' | awk -v a=$s -v b=$3 '{print b","a",P50,"$1}'
    awk -v a=$s -v b=$3 '/completeness/ {print b","a",completeness,"$2}' $1_${s}_$2.merqury.full.stats
    awk -v a=$s -v b=$3 '/QV/ {print b","a",QV,"$2}' $1_${s}_$2.merqury.full.stats
    awk -v a=$s -v b=$3 'NR<31 {c+=$2} END {print b","a",gaps,"c}' $1_${s}_$2.gaps.txt
    awk -v a=$s -v b=$3 '/Complete BUSCOs/ {print b","a",BUSCO,"$1}' $1_${s}_$2.BUSCO.txt
done
