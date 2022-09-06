##########################################
#
#       A.altamirani skin microbiome          
#       Martínez-Ugalde et al. 2022
#           emartug@gmail.com
#
##########################################


#Core features for metamorphic and non-metamorphic samples in QIIME2

qiime feature-table filter-samples \
--i-table Altamirani_core-metrics-results/rarefied_table.qza \
--m-metadata-file metadata_altamirani.txt \
--p-where "SampleType IN ('Adult')" \
--o-filtered-table met90-table.qza

qiime feature-table summarize \
--i-table met90-table.qza \
--o-visualization met90-table.qzv

qiime feature-table core-features \
--i-table met95-table.qza \
--p-min-fraction 0.90 \
--p-max-fraction 1 \
--p-steps 10 \
--o-visualization met90-core-table.qzv

qiime feature-table filter-samples \
--i-table Altamirani_core-metrics-results/rarefied_table.qza \
--m-metadata-file metadata_altamirani.txt \
--p-where "SampleType IN ('Juvenile')" \
--o-filtered-table nmet90-table.qza

qiime feature-table summarize \
--i-table nmet90-table.qza \
--o-visualization nmet90-table.qzv

qiime feature-table core-features \
--i-table  nmet90-table.qza \
--p-min-fraction 0.90 \
--p-max-fraction 1 \
--p-steps 10 \
--o-visualization nmet90-core-table-rare.qzv

qiime feature-table filter-samples \
--i-table Altamirani_core-metrics-results/rarefied_table.qza \
--m-metadata-file metadata_altamirani.txt \
--p-where "SampleType IN ('Sediment')" \
--o-filtered-table sed90-table.qza

qiime feature-table summarize \
--i-table sed90-table.qza \
--o-visualization sed90-table.qzv

qiime feature-table core-features \
--i-table  sed-table.qza \
--p-min-fraction 0.90 \
--p-max-fraction 1 \
--p-steps 10 \
--o-visualization sed90-core-table-rare.qzv

qiime feature-table filter-samples \
--i-table Altamirani_core-metrics-results/rarefied_table.qza \
--m-metadata-file metadata_altamirani.txt \
--p-where "SampleType IN ('Water')" \
--o-filtered-table wat90-table.qza

qiime feature-table summarize \
--i-table wat90-table.qza \
--o-visualization wat90-table.qzv

qiime feature-table core-features \
--i-table  wat90-table.qza \
--p-min-fraction 0.90 \
--p-max-fraction 1 \
--p-steps 10 \
--o-visualization wat90-core-table-rare.qzv
