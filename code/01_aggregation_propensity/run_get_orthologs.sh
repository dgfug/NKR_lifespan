#!/bin/bash
# Author: Savandara Besse
# Created: 03-15-2018

##################################

#### Command line to generate inparanoid outputs
# perl inparanoid /media/DATA1/savvy/inparanoid/HG-MM_Uniprot/MM /media/DATA1/savvy/inparanoid/HG-MM_Uniprot/HG

### Custom orthologs with last proteome version - Last run 04-11-2019
## Parameters
sqltable='/media/DATA1/savvy/inparanoid/HG-MM_Uniprot/sqltable.MM-HG'
output='_orthologs.faa'

### For h_glaber
proteome='/media/DATA1/savvy/BLAST/model/het_glab/het_glab'
filter='HG'
org='HG'


# ### For m_musculus
# proteome='/media/DATA1/savvy/BLAST/model/mus_musc/mus_musc'
# filter='MM'
# org='MM'

/usr/bin/python3 in_getOrthologs.py -x $filter -p $proteome -t $sqltable -o ./uniprot_results/$org$output &
