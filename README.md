# Inserted genomic element in mammalian genomes

1) I'll scan liftOver chain files (hg38 based) to find insertions in other genomes relative to human

2) insertions scan output for a given assembly is a bed file with four fields, eg:
	in the file: output_examples/mm10.gapsInChains.bed 
```
	chr5	178604087	178604087	196_SS_gap587_0_13
	chr7	154699597	154699710	152_DS_gap3018_113_99
	chr13	78498645	78498749	23_DS_gap22218_104_1451
	chr20	9824436	9824623	22_DS_gap8377_187_119
```	
   the 4'th field string read as:	chainI_gap-type_gap-number_ref.-deleted_query_insertions

   chain_ID		// the ID of the liftOver chain (downloaded from UCSC genome browser) on which we found a gap when comparing the human genome (hg38) to a query genome

   gap_type		// SS if there is insertion in the query genome and no deletion in the reference; DS is we have ref deletion and query insertion (i.e., Single - and doble sided gaps, respectively)

   <ref.deleted>	// size of chain-gap in the reference genome (i.e., gaps on the target species chain)

   <query insertions>	// size of inserted bases in the target species
