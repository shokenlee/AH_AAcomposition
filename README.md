# Amphipahtic helices_Amino-acid composition

Scripts for amino acid composition analysis of amphipathic helices in [Lee et al. bioRxiv 2024](https://www.biorxiv.org/content/10.1101/2024.11.14.623600v2)

Tested on Python v3.8 and v3.11 with Anaconda v2.6.

Amino acid composition analysis in Figure 2 and NEES analysis in Figure 5 are done in paralell.

## Files
### Source
- `NEES.xlsx`: NE enrichment sccore data from Figure 4
- `AH_sequences.xlsx`: AH names and amino acid sequence
- `ddF.csv`: ddF values
### Analysis
- `NEES_Import.ipynb`: NE enrichment scores were analyzed to give their medians
  - Output: `NEES_summary.csv` and `NEES_Normalized.csv`
- `NEESvsAH_composition.ipynb`: Amino acid composition analysis in hydrophobic and hydrophilic faces and the results were merged with NE enrichment score medians
  - Output: `NEES_summary_with_heliquest.csv`
- `NEES_Transform.ipynb`: Box-Cox transformation of NE enrichment score medians
  - Output: `NEES_summary_with_heliquest_Median_Transformed.csv`
- `Heatmap.ipynb`: Heatmap anlaysis
  - Output: `All.pdf` and `heatmap.csv` 
- `LDA.ipynb`: Linear discriminant analysis
  - Output: `lda_plot.pdf` and `lda_coefficients.csv` and `lda_results.csv`
- `ddF_plot.ipynb`: Ploting NE enrichment score with ddF values
  - Output: `ddFvsNEES.csv`
### Module
- `module_heliquest_like.py`: functions and class for parameter disection of protein helix

