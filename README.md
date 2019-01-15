# ptm_hotspot: prediction of PTM hotspots

## Description

Method to infer PTM (currently only phosphorylation) hotspots in conserved protein alignments based on (Strumillo *et al*. *bioRxiv*).

## Dependencies

`ptm_hotspot` requires Python (v3) as well as the following packages: 

- numpy
- pandas
- scipy
- statsmodels

## How to predict hotspots

If you want to predict all hotspots in all protein domains just run:

`python3 ptm_hotspots.py -o predicted_hotspots.csv`

To obtain particular domain predictions:

`python3 ptm_hotspots.py -o kinase_domain_hotspots.csv -d PF00069`

To obtain hotspot residue predictions instead of hotspot ranges 

`python3 ptm_hotspots.py -o hotspot_residues.csv --printSitePredictions`

You can obtain more help and options by typing:

`python3 ptm_hotspots.py -h`

Note: Since the Bonferroni correction depends on the total number of predictions, small disimilarities might emerge in the same domain hotspots depending on whether you run only a domain or the full set of domains. Similarly, the stochastic nature of the permutation analysis might make the results vary between runs.

### Customizing alignments or PTM data

By default `ptm_hotspot` uses a database containing precalculated domain alignments (as described in Strumillo et al.) as well as a collection of phosphorylated residues derived from public high-throughput mass spectrometry experiments. In order to update the database please consider the next points:

##### Alignments

Every alignment file should be in FASTA format and the header should contain the start and the end of the domains in the alignment coordinates separated by ";". For example:
	
	>EDP05298 pep:known supercontig:v3.1:DS4 ;51;337

For full protein predictions just include the first and last positions in the multiple sequence alignment.

##### PTM database

The ptm database should be included as a `csv` file containing id, amino acid and position of the phosphosite within the protein.

## Citation

- Strumillo, M. J., Oplova, M., Viéitez, C., Ochoa, D., Shahraz, M., Busby, B. P., et al. (2018). Conserved phosphorylation hotspots in eukaryotic protein domain families. bioRxiv. [https://doi.org/10.1101/391185](https://doi.org/10.1101/391185)
