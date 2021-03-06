# FAECES*: Fast and Accurate Enrichment of Canid Excrement for Species (\*and other analyses)  

This repository contains scripts to analyze noninvasive DNA from Mojave desert canids using the FAECES* assay.  

## Citation  
Please cite:  
Parker, L.D., Campana, M.G., Quinta, J.D., Cypher, B., Rivera, I., Fleischer, R.C., Ralls, K., Wilbert, T.R., Boarman, R., Boarman, W.I. & Maldonado, J.E. 2022. An efficient method for simultaneous species, individual, and sex identification via in-solution SNP capture of low-quality scat samples. *Molecular Ecology Resources*. 22: 1345-1361. DOI: [10.1111/1755-0998.13552](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13552).  

## Prerequistes  
`canid_sex.rb` requires [Ruby](www.ruby-lang.org) >= 2.4.  

## Installation  
Clone the repository: `git clone https://github.com/campanam/FAECES`  
Make the scripts executable: `chmod +x *.rb`  

## Usage  
`canid_sex.rb` performs sex determination using a previously called SNP VCF of the *ZFX/ZFY* region. Options include:  

    -i, --input [FILE]       Input VCF file of SNPs in the ZFX/ZFY region. *Required*.
    -s, --species [FILE]     Species assignment file. *Optional and Experimental*.
    -k, --males [FILE]       List of known males to empirically calculate Y-allele frequency. *Optional*.
        --mean               Infer Y frequency using mean Y frequency across known male samples
    -z, --zscore [VALUE]     Infer mean Y frequency from known male samples and set z-score bound
        --min                Infer Y frequency using minimum Y frequency across known male samples
    -m, --minalleles [VALUE] Minimum number of allele reads to call sex (Default = 10)
    -y, --yfreq [VALUE]      Expected Y-allele frequency (Default = 0.5). *Required if not inferred from known males*.
    -a, --alpha [VALUE]      Alpha value for statistical significance (Default = 0.05).
        --yx [VALUE]         Use specified Y-X allele ratio to determine sex rather than statistical test.
    -Z [VALUE]               Use z-score cut-off to determine Y-X ratio.
    -h, --help               Show help

The males file is a simple text list. See [canid_males_ZFX.txt](canid_sex_data/canid_males_ZFX.txt) for an example.  
The species file is a CSV in the format `sample,species`. Use 'KF' to denote kit fox and use kit-fox-specific Y alleles. *NB: this option is experimental and has not been fully tested or debugged.*  

The final sex determination datasets are available [here](canid_sex_data). To obtain the raw sex determinations, use the following commands:  

Coyotes: `ruby canid_sex.rb -i Coy_ZFX.raw.vcf -k canid_males_ZFX.txt --min > Coy_ZFX_sex.tsv`  
Kit foxes: `ruby canid_sex.rb -i KitFox_ZFX.raw.vcf -k canid_males_ZFX.txt --min > KitFox_ZFX_sex.tsv`  


