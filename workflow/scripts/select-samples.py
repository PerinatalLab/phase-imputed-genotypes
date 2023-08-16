import pandas as pd
import numpy as np

#### Read data

d= pd.read_csv(snakemake.input[0], sep= ' ', header= 0, names= ['FID', 'IID', 'dad', 'mom', 'sex', 'pheno'])

#### Reformat kin file so that all individuals appear in ID1 column

kin= pd.read_csv(snakemake.input[1], header= 0, sep= '\t')

kin1= kin.copy()

kin1.columns= ['FID', 'ID2', 'ID1', 'N_SNP', 'Z0', 'Phi', 'HetHet', 'IBS0', 'HetConc', 'HomIBS0', 'Kinship', 'IBD1Seg', 'IBD2Seg', 'PropIBD', 'InfType', 'Error']

kin= pd.concat([kin, kin1])
kin= kin.loc[kin.InfType== 'PO', :]

#### Parent-offspring relationships  // Trios

trios= d.loc[(d.dad.str.contains('_')) & (d.mom.str.contains('_')) & (d.mom.isin(d.IID)) & (d.dad.isin(d.IID)), :]

kin_trios= kin.loc[kin.ID1.isin(trios.IID), :]

# Check that it is actually parent-offspring trios

temp_trios= pd.merge(trios, kin_trios, left_on= 'IID', right_on= 'ID1')
temp_trios_dad= temp_trios.loc[temp_trios.ID2== temp_trios.dad, :]
temp_trios_mom= temp_trios.loc[temp_trios.ID2== temp_trios.mom, :]
trios= trios.loc[(trios.IID.isin(temp_trios_mom.IID.values)) & (trios.IID.isin(temp_trios_dad.IID.values)), :]

# Duos

dad_duos= d.loc[(d.dad.str.contains('_')) & (~d.mom.str.contains('_')), :]
mom_duos= d.loc[(~d.dad.str.contains('_')) & (d.mom.str.contains('_')), :]
kin_dad_duos= pd.merge(kin, dad_duos, left_on= ['ID1', 'ID2'], right_on= ['IID', 'dad'])
kin_mom_duos= pd.merge(kin, mom_duos, left_on= ['ID1', 'ID2'], right_on= ['IID', 'mom'])

dad_duos= dad_duos.loc[dad_duos.IID.isin(kin_dad_duos.IID.values), :]
mom_duos= mom_duos.loc[mom_duos.IID.isin(kin_mom_duos.IID.values), :]
pedigree= pd.concat([trios, dad_duos, mom_duos])

related_ids= pd.concat([pedigree[['FID', 'IID']], pedigree[['FID', 'dad']], pedigree[['FID', 'mom']]])
related_ids['IID']= np.where(related_ids.IID.isna(), related_ids.dad, related_ids.IID)
related_ids['IID']= np.where(related_ids.IID.isna(), related_ids.mom, related_ids.IID)
related_ids= related_ids.loc[related_ids.IID.str.contains('_'), :]
related_ids.drop_duplicates(['FID', 'IID'], inplace= True, keep= 'first')

pedigree['dad']= np.where(~pedigree.dad.str.contains('_'), 'NA', pedigree.FID + '_' + pedigree.dad)
pedigree['mom']= np.where(~pedigree.mom.str.contains('_'), 'NA', pedigree.FID + '_' + pedigree.mom)
pedigree['IID']= pedigree.FID + '_' + pedigree.IID

# non related individuals - aka, individuals without parents/ offspring

sind= d.loc[~d.IID.isin(related_ids.IID.values), :]

sind.to_csv(snakemake.output[1], sep= '\t', header= False, index= False, columns= ['FID', 'IID'])
related_ids.to_csv(snakemake.output[0], sep= '\t', header= False, index= False, columns= ['FID', 'IID'])


related_ids['IID']= related_ids.FID + '_' + related_ids.IID
sind['IID']= sind.FID + '_' + sind.IID

pedigree.to_csv(snakemake.output[2], sep= '\t', header= False, index= False, columns= ['IID', 'dad', 'mom'])
