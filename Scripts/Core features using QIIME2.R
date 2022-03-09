##########################################
#
#       A.altamirani skin microbiome          
#       Martínez-Ugalde et al. 2022
#           emartug@gmail.com
#
##########################################


#Core features for metamorphic and non-metamorphic samples in QIIME2
qiime feature-table filter-samples \
--i-table filtrar-table.qza \
--m-metadata-file metadata_altamirani.txt \
--p-where "SampleType IN ('Metamorphic')" \
--o-filtered-table met-table.qza

qiime feature-table summarize \
--i-table adult-table.qza \
--o-visualization met-table.qzv 


qiime feature-table core-features \
--i-table met-table.qza \
--p-min-fraction 0.8 \
--p-max-fraction 1 \
--p-steps 20 \
--o-visualization met-core-table.qzv

qiime feature-table filter-samples \
--i-table filtrar-table.qza \
--m-metadata-file metadata_altamirani.txt \
--p-where "SampleType IN ('Non-metamorphic')" \
--o-filtered-table nmet-table.qza

qiime feature-table summarize \
--i-table nmet-table.qza \
--o-visualization nmet-table.qzv

qiime feature-table core-features \
--i-table  juvenile-table.qza \
--p-min-fraction 0.8 \
--p-max-fraction 1 \
--p-steps 20 \
--o-visualization nmet-core-table-rare.qzv