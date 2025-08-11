# Bidirectional global source-sink dynamics of dengue viruses and epidemic establishment in non-endemic areas

## Before starting
1. All code and data contained within this repository are released under the CC BY-NC-SA License. 
2. Because of GISAID [terms of use](https://www.gisaid.org/registration/terms-of-use/), sequence data are not provided. Sequences can be downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) or [GISAID](https://www.gisaid.org/).

## Repository structure and usage
The structure of this repository is shown below. 
```
├── 1. Global analysis (Fig 2)
│   ├── phylogenetic_analysis
│   └── phylogeographic_analysis
├── 2. Annual analysis (Fig 3)
│   ├── phylogenetic_analysis
│   └── phylogeographic_analysis
├── 3. China analysis (Fig 4)
│   ├── phylogenetic_analysis
│   └── phylogeographic_analysis
├── 4. GLM analysis (Fig 5)
│   ├── phylogenetic_analysis
│   └── phylogeographic_analysis
├── Acknowledgement_table
├── README.md
```

# Phylogenetic analyses
Genetic sequences have been scrubbed from BEAST XML files because of GISAID [terms of use](https://www.gisaid.org/registration/terms-of-use/).

# Phylogeographic analyses
Before running the XML files, a file of empirical posterior trees needs to be obtained and the filenames of the empirical tree files in the XML files may need to be modified. Such an empirical tree file can be obtained from the analyses inside `phylogenetic_analysis/` folder.

<h1> License </h1>
<h4>CC BY-NC-SA 4.0 </h4>

Attribution-NonCommercial-ShareAlike 4.0 International
This license requires that reusers give credit to the creator. It allows reusers to distribute, remix, adapt, and build upon the material in any medium or format, for noncommercial purposes only. If others modify or adapt the material, they must license the modified material under identical terms.

BY: Credit must be given to you, the creator.

NC: Only noncommercial use of your work is permitted.
Noncommercial means not primarily intended for or directed towards commercial advantage or monetary compensation.

SA: Adaptations must be shared under the same terms.

Copyright (c) 2025 Zhiyuan Chen
