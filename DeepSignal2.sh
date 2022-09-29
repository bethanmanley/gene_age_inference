
##################################################################
#Concatenate *.fastq.gz files from all runs
##################################################################

cat 20210412_1434_MN32650_FAP96969_01caaf79/guppy_out/*.fastq.gz 20201029_0937_MN32650_FAO25988_012182a7/guppy_out/*.fastq.gz 20201007_1629_MN32650_FAO66202_c64dbc4b/guppy_out/*.fastq.gz > ~/mfre/guppy_out_pool/guppy_pool.fastq.gz
wc -c guppy_pool.fastq.gz
# Rhirr.DAOM197198: 12,002,780,365 bytes
gunzip -c guppy_pool.fastq.gz > guppy_pool.fastq
wc -c guppy_pool.fastq
# Rhirr.DAOM197198: 25,458,626,525 bytes


##################################################################
# Pool all fast5s into newdir and convert multi fast5s to single fast5s
##################################################################

# Pool all multi fast5s together into /mfre/multi_fast5_pool (include both pass and fail reads)
mv 20210412_1434_MN32650_FAP96969_01caaf79/fast5/*.fast5 ~/multi_fast5_pool
mv 20201029_0937_MN32650_FAO25988_012182a7/fast5/*.fast5 ~/multi_fast5_pool
mv 20201007_1629_MN32650_FAO66202_c64dbc4b/fast5/*.fast5 ~/multi_fast5_pool

# Convert multi to single
multi_to_single_fast5 -i multi_fast5_pool -s single_fast5_pool -t 30 --recursive


##################################################################
# tombo preprocess
##################################################################

tombo preprocess annotate_raw_with_fastqs --fast5-basedir single_fast5_pool --fastq-filenames guppy_out_pool/guppy_pool.fastq --basecall-group Basecall_1D_000 --basecall-subgroup BaseCalled_template --overwrite --processes 10

##################################################################
# tombo resquiggle
##################################################################

# Working with the whole assembly, including mito.genome and possibly endobacterial.symbiont/contamination
export HDF5_PLUGIN_PATH=/usr/local/hdf5/lib/plugin
tombo resquiggle single_fast5_pool Rirr.curated_primary.no_mt.unscrubbed.fa --processes 10 --corrected-group RawGenomeCorrected_000 --basecall-group Basecall_1D_000 --overwrite --ignore-read-locks


##################################################################
# deepsignal2 extract features
##################################################################

/usr/local/bin/deepsignal2 extract -i single_fast5_pool \
--reference_path Rirr.curated_primary.no_mt.unscrubbed.fa \
-o deepsignal2_out/Rirr.manley.fast5s.CG.features.tsv \
--corrected_group RawGenomeCorrected_000 \
--nproc 15 \
--motifs CG


##################################################################
# Call modifications
##################################################################
CUDA_VISIBLE_DEVICES=0 deepsignal2 call_mods --input_path single_fast5_pool \
--model_path ~/mfre/models/model.dp2.CG.R9.4_1D.human_hx1.bn17_sn16.both_bilstm.b17_s16_epoch4.ckpt \
--result_file ~/deepsignal2_out/Rirr.manley.fast5s.CG.call_mods.tsv \
--corrected_group RawGenomeCorrected_000 \
--reference_path Rirr.curated_primary.no_mt.unscrubbed.fa \
--motifs CG \
--nproc 10 \
--nproc_gpu 3


##################################################################
# Modification frequency output file
##################################################################
# I copied the call_modification_frequency.py and txt_formater.py scripts from github and copied them in my path (~/nanopore/)

python call_modification_frequency.py --input_path deepsignal2_out/Lyc_spec.CG.call_mods.tsv \
--result_file deepsignal2_out/Rhizophagus_irregularis_DAOM197198_mCG_mods_frequency.tsv \
--prob_cf 0

#Sort the output. (I deleted the non-sorted one and removed 'sorted' from the sorted one)
#sort -V -k1,1 -k2,2 *.tsv > *_sorted.tsv
