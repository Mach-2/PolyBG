import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
np.random.seed(0)
# Read in data from GWAS Catalog 
traits = pd.read_csv('traits-clean.csv')

# Read in GIANT GWAS 
giant_gwas = pd.read_table('GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.tsv')[["SNPID","EFFECT_ALLELE_FREQ","P"]]
giant_PGS = pd.read_table('GIANT_HEIGHT_YENGO_2022_PGS_WEIGHTS_ALL.tsv')[["SNPID","PGS_EFFECT_ALLELE","PGS_WEIGHT"]]
joined = giant_gwas.set_index('SNPID').join(giant_PGS.set_index('SNPID'),how='inner')
# print(joined.head())
# print("Number of traits: ",len(joined))

filtered = joined[joined['P'] < 0.00000005] # GIANT study said genome-wide significance defined here as P<5x10^-8
# print("Traits with P < 5 x 10^-8: ", len(filtered))

# print("Sorted by effect size: ", filtered.sort_values('PGS_WEIGHT').tail(), filtered.sort_values('PGS_WEIGHT').head())

traits = filtered.rename(columns={'EFFECT_ALLELE_FREQ': 'RAF', 'PGS_WEIGHT': 'Effect'})


# effects = list(traits['Effect'])
# effects.sort()
# # plt.plot(effects)
# sns.distplot(effects)
# # plt.axhline(0.5)
# plt.show()