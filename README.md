# SC2 Mutation Frequency Calculator (SMFC)
A covSonar Utility tool to detect characteristic and signature mutations based on mutation profiles. 
This tool allows the calculation of mutation frequencies for a user-defined timeframe for specific lineages. 

## 1. Setup 

### 1.1 Prerequisites
Conda and python are necessary.
[CovSonar](https://github.com/rki-mf1/covsonar) is the key tool to extract mutation profiles from. It has a database with consensussequences sequenced in the course of the pandemic assigned to lineages which can be queried. 

### 1.2 Install
Proceed as follows to install SMFC:
```bash
# download the repository to the current working directory using git 
git clone https://github.com/rki-mf1/sc2-mutation-frequency-calculator.git
# build the custom software environment using conda [recommended]
conda env create -n sonar -f sc2-mutation-frequency-calculator/smfc.env.yml
# activate the conda evironment if built 
conda activate smfc
```

### 1.3 Options/--help

| option              | value(s)                                                              | note |
|---------------------|-----------------------------------------------------------------------|------| 
| --acc               | one or more genome accessions (e.g. NC_045512.2)                      |      |
| --lineage           | one or more pangolin lineages (e.g. B.1.1.7)                          |      |
| --zip               | one or more zip codes (e.g. 10627)                                    | zip codes are dynamically extended to the right side, e.g. 033 matches to all zip codes starting with 033|
| --date              | one or more dates or date ranges (e.g. 2021-01-01)                    | single dates are formatted as YYYY-MM-DD while date ranges can be defined by YYYY-MM-DD:YY-MM-DD (from:to) |
| --submission_date   | one or more dates or date ranges (e.g. 2021-01-01)                    | single dates are formatted as YYYY-MM-DD while date ranges can be defined by YYYY-MM-DD:YY-MM-DD (from:to) |
| --lab               | one or more labs (e.g. L1)                                            |      |
| --source            | one or more data sources (e.g. DESH)                                  |      |
| --collection        | one or more data collections (e.g. RANDOM)                            |      |
| --technology        | one or more sequencing technologies (e.g. Illumina)                   |      |
| --platform          | one or more sequencing platforms (e.g. MiSeq)                         |      |
| --chemistry         | one or more sequencing chemistries (e.g. Cleanplex)                   |      |
| --software          | one software tool used for genome reconstruction (e.g. covPipe)       |      |
| --version           | one software tool version used for genome reconstruction (e.g. 3.0.5) | needs --software defined |
| --material          | one or more sample materials (e.g. 'nasal swap')                      |      |
| --min_ct            | minimal ct value (e.g. 20)                                            |      |
| --max_ct            | maximal ct value (e.g. 20)                                            |      |
 

### 1.4 Default
Mutation frequency matrix 
(figure)

*Parent Lineage: 
*Number of sequences detected:
*Labdiversity:

How a frequency matrix can be created: 
```bash
# download the repository to the current working directory using git 
git clone https://github.com/rki-mf1/sc2-mutation-frequency-calculator.git
# build the custom software environment using conda [recommended]
conda env create -n sonar -f sc2-mutation-frequency-calculator/smfc.env.yml
# activate the conda evironment if built 
conda activate smfc
```

## 2. Examples

### 2.1 variant-specific pcr test-design
signature mutations can be used to determine which mutations acuretly define a lineage (nothing else) which can be used to test for specific lineage in a variant-specific pcr test-design:
(figure)

How signature mutations can be calculated:  

```bash
# download the repository to the current working directory using git 
git clone https://github.com/rki-mf1/sc2-mutation-frequency-calculator.git
# build the custom software environment using conda [recommended]
conda env create -n sonar -f sc2-mutation-frequency-calculator/smfc.env.yml
# activate the conda evironment if built 
conda activate smfc
```

### 2.2 consensus^2 for representing a lineage

How signature mutations can be calculated:  

```bash
# download the repository to the current working directory using git 
git clone https://github.com/rki-mf1/sc2-mutation-frequency-calculator.git
# build the custom software environment using conda [recommended]
conda env create -n sonar -f sc2-mutation-frequency-calculator/smfc.env.yml
# activate the conda evironment if built 
conda activate smfc
```

## 3. Best practice

### 3.1 Inputformat

| Lineage  |  | Mutation  | 
| ---------|  | ---------|
| BA.5.2.1 |  | S:D46Y |
| BA.4.6   |  | ORF1ab:S367G   |
| BE.1     |  | N:A679D   |

### 3.2 Presence of overlap in the data 
See https://github.com/rki-mf1/sc2-mutation-frequency-calculator/issues/4#issue-1755232495

## 4. Contribution 
covSonar has been very carefully programmed and tested, but is still in an early stage of development. You can contribute to this project by reporting problems or writing feature requests to the issue section under https://github.com/rki-mf1/sc2-mutation-frequency-calculator/-/issues

Your feedback is very welcome!

## 5. FAQ

### characteristic vs signature mutations? 
figure

### consensus^2?
Other then the conventional consensus method consensus^2 is used to build a robust and represenative consensus of a number of samples for a lineage by introducing the most frequent mutations (default cut-off:10) in a timeframe in a reference genome (default: Wuhan). 

Could be used for primer selection (and building phylogenetic trees):
figure






