# PALs_tn5seq_VariantCall

# align reads to  MSH2_1-934_NNN.fa using bbmap

export BD=/nfs/kitzman2/isaac/miseq/20180419_miseq/
export DATE=`date "+%Y%m%d"`
cd $BD
mkdir $BD/${DATE}MSH2_1-934_NNN_out
export OUT=$BD/${DATE}MSH2_1-934_NNN_out

export datadir=/nfs/kitzman2/isaac/miseq/20180419_miseq/clip_fqt/

export pfx=${datadir}XJ_XJ
is=(${pfx}*.read1.fq.gz)
js=(${pfx}*.read3.fq.gz)

for ((i = 0; i < 5; i++)); do /nfs/kitzman2/isaac/softwares/bbmap/bbmap.sh in1=${is[i]} in2=${js[i]} out=$OUT/${is[i]#$datadir}.sam ref=/nfs/kitzman2/isaac/miseq/refs/MSH2_1-934_NNN.fa nodisk -Xmx30g semiperfectmode=t; done


cd $OUT
export suffix=.read1.fq.gz.sam
for i in *.sam; do samtools view -bS $i | samtools sort > ${i%$suffix}.bam; echo ${i%$suffix}_completed; done

rm *.sam
for i in *.bam; do samtools index $i; done



# Call variants using the 20180117_VarCall_todf.py (this could take half an hour)

#please change the first arguement to your actual path, eg  bd=/nfs/kitzman2/isaac/miseq/20180419_miseq/20180424MSH2_1-934_NNN_out/
for i in *.bam; do python /nfs/kitzman2/isaac/xj_scripts/20180412_VarCall_todf_plasmidref.py /nfs/kitzman2/isaac/miseq/20180419_miseq/20180424MSH2_1-934_NNN_out/ $i &
done

