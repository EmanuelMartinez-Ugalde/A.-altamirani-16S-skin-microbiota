##########################################
#
#       A.altamirani skin microbiome          
#       Martínez-Ugalde et al. 2022
#           emartug@gmail.com
#
##########################################

#ANCOM between metamoprhic and non-metamophic samples at family level

qiime feature-table filter-samples \
--i-table Altamirani_filtered_OTU-table.qza \
--m-metadata-file metadata_altamirani.txt \
--p-where "SampleType IN ('Metamorphic','Non-metamorphic')" \
--o-filtered-table ANCOM/met-nmet-table.qza

qiime feature-table filter-features \
--i-table  ANCOM/met-nmet-juvenile-table.qza \
--p-min-frequency 50 \
--o-filtered-table ANCOM/met-nmet-feature-frequency-filtered-table.qza

qiime taxa collapse \
--i-table ANCOM/met-nmet-feature-frequency-filtered-table.qza \
--i-taxonomy  Altamirani_taxonomy.qza \
--p-level 5 \
--o-collapsed-table ANCOM/met-nmet-table-coll.qza

qiime composition add-pseudocount \
--i-table ANCOM/met-nmet-table-coll.qza \
--o-composition-table ANCOM/comp-met-nmet-table.qza

qiime composition ancom \
--i-table  ANCOM/comp-met-nmet-table.qza \
--m-metadata-file metadata_altamirani.txt \
--m-metadata-column SampleType \
--o-visualization ANCOM/L5-met-nme.qzv

#ANCOM for matamorph9ic samples between winter-spring

qiime feature-table filter-samples \
--i-table Altamirani_filtered_OTU-table.qza \
--m-metadata-file metadata_altamirani.txt \
--p-where "SampleType IN ('Metamorphci')" \
--o-filtered-table ANCOM/met-table.qza

qiime feature-table filter-samples \
--i-table ANCOM/met-table.qza \
--m-metadata-file metadata_altamirani.txt \
--p-where "Season IN ('Winter','Spring')" \
--o-filtered-table ANCOM/met-winter-spring-table.qza

qiime feature-table filter-features \
--i-table ANCOM/met-winter-spring-table.qza \
--p-min-frequency 50 \
--o-filtered-table ANCOM/met-winter-spring-feature-frequency-filtered-table.qza

qiime taxa collapse \
--i-table  ANCOM/met-winter-spring-feature-frequency-filtered-table.qza \
--i-taxonomy Altamirani_taxonomy.qza \
--p-level 5 \
--o-collapsed-table ANCOM/met-winter-spring-table.qza

qiime composition add-pseudocount \
--i-table ANCOM/met-winter-spring-table.qza \
--o-composition-table ANCOM/comp-met-winter-spring-table.qza

qiime composition ancom \
--i-table ANCOM/comp-met-winter-spring-table.qza \
--m-metadata-file metadata_altamirani.txt \
--m-metadata-column Season \
--o-visualization ANCOM/L5-met-winter-spring.qzv

#ANCOM for non-metamorphic samples between autummn-winter and winter-spring

qiime feature-table filter-samples \
--i-table Altamirani_filtered_OTU-table.qza \
--m-metadata-file metadata_altamirani.txt \
--p-where "SampleType IN ('Non-metamorphic')" \
--o-filtered-table ANCOM/nmet-table.qza

qiime feature-table filter-samples \
--i-table ANCOM/nmet-table.qza \
--m-metadata-file metadata_altamirani.txt \
--p-where "Season IN ('Autumn','Winter')" \
--o-filtered-table ANCOM/nmet-autmn-winter-table.qza

qiime feature-table filter-features \
--i-table ANCOM/nmet-autmn-winter-table.qza \
--p-min-frequency 50 \
--o-filtered-table ANCOM/nmet-autumn-winter-feature-frequency-filtered-table.qza

qiime taxa collapse \
--i-table ANCOM/nmet-autumn-winter-feature-frequency-filtered-table.qza \
--i-taxonomy  Altamirani_taxonomy.qza \
--p-level 5 \
--o-collapsed-table ANCOM/juvenile-autumn-winter-table.qza

qiime composition add-pseudocount \
--i-table ANCOM/nmet-autumn-winter-table.qza \
--o-composition-table ANCOM/comp-nmet-autumn-winter-table.qza

qiime composition ancom \
--i-table ANCOM/comp-nmet-autumn-winter-table.qza \
--m-metadata-file Metadata/metadata_altamirani.txt \
--m-metadata-column Season \
--o-visualization ANCOM/L5-ancom-nmet-autumn-winter.qzv
##############################################################################################################
qiime feature-table filter-samples \
--i-table ANCOM/nmet-table.qza \
--m-metadata-file metadata_altamirani.txt \
--p-where "Season IN ('Winter','Spring')" \
--o-filtered-table ANCOM/juvenile-winter-spring-table.qza

qiime feature-table filter-features \
--i-table ANCOM/nmet-winter-spring-table.qza \
--p-min-frequency 50 \
--o-filtered-table ANCOM/nmet-winter-spring-feature-frequency-filtered-table.qza

qiime taxa collapse \
--i-table ANCOM/nmet-winter-spring-feature-frequency-filtered-table.qza \
--i-taxonomy  Altamirani_taxonomy.qza \
--p-level 5 \
--o-collapsed-table ANCOM/nmet-winter-spring-table-L4.qza

qiime composition add-pseudocount \
--i-table ANCOM/nmet-winter-spring-table-L4.qza \
--o-composition-table ANCOM/comp-nmet-winter-spring-table-L4.qza

qiime composition ancom \
--i-table ANCOM/comp-nmet-winter-spring-table-L4.qza \
--m-metadata-file Metadata/metadata_altamirani.txt \
--m-metadata-column Season \
--o-visualization ANCOM/L5-ancom-nmet-winter-spring.qzv