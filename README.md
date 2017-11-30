# LMSACbacteria
In this work, we developed a bacterial MALDI TOF spectra identification approach without using standard spectra library, but only protein sequences from Uniprot database. Genetic algorithm with double cross-validation was performed to get the proteins encoded by the most informative genes as protein panels to be considered for bacterial identification. At genus level, the identification accuracy exceeded 80%. Although the performance is still behind to the performance of library search based approach, it can be used as a complementary method when standard library is not available. It can also be used in diverse application areas such as clinical microbiology, biodefense, food safety and environmental health.
## Requirement
Python 3.5 with `numpy`, `pandas`, and `matplotlib` 

[Anaconda] (https://anaconda.org/anaconda/python) is recommended.

## Database
### Spectra
`LMSACbacteria` uses peak lists for bacteria identification. A peak list is stored in a space-separated text (.txt) file including the m/z and intensity of all the peaks from a mass spectrum of a strain of bacteria with the first column as m/z and the second column as intensity, e.g.:
```
4025.85 14.9166578754281
4179.39 10.0793564012325
4195.16 100
4279.94 13.6100699949239
4296.99 13.759296917905
4339.2 16.2279091177302
...
9998.94 4.53344723468826
10103.51 25.6021474910836
10245.55 14.5594503134826
10457.14 9.18592847464415
10929.3 2.39529317439903
11115.92 10.3301379804275
11339.22 8.75138710039803
11702.38 2.53578789353994
```
- European Consortium of Microbial Resources Centres (EM-baRC): [http://www.embarc.eu/deliverables/EMbaRC_D.JRA2.1.4_D15.33_MALDI-TOF_DB.pdf]

- Robert Koch-Institute: [http://www.microbe-ms.com/]

- The Public Health Agency of Sweden: [http://spectra.folkhalsomyndigheten.se/spectra/]
### Uniprot
`LMSACbacteria` uses proteomic sequences in **Uniprot** to match the MALDI-TOF spectra.

```
>tr|A0A087INJ0|A0A087INJ0_VIBVL UDP-N-acetylmuramoyl-L-alanyl-D-glutamate--2,6-diaminopimelate ligase OS=Vibrio vulnificus GN=murE PE=3 SV=1
MRNTMNLTNLLAPWLDCPELADITVQSLELDSRQVKQGDTFVAIVGHVVDGRQYIEKAIE
QGANAIIAQSCQQYPSGLVRYQQNVVIVYLEKLDEKLSQLAGRLYQHPEMSLIGVTGTNG
KTTITQLIAQWLELAGQKAAVMGTTGNGFLNALQPAANTTGNAVEIQKTLADLQQQGAKA
TALEVSSHGLVQGRVKALQFAAGVFTNLSRDHLDYHGTMEAYAQAKMTLFTEHQCQHAII
NLDDEVGAQWFQELKQGVGVSLYPQDASVKALWASSVAYAESGITIEFEGCFGQGRLHAP
LIGEFNATNLLLALATLLALGVDKQALLGSAANLRPVLGRMELFQVNSKAKVVVDYAHTP
DALEKALQALRVHCTGHLWAIFGCGGDRDKGKRPMMAEIAERLADHVVLTDDNPRSEDPA
MIVQDMLAGLTRADSAVVEHDRFSALQYALDNAQADDIILLAGKGHEDYQVLKHQTVHYS
DRESAQQLLGISS
```
## Usage
Run [main.py] to load the functions.
### Identificaiton

### Cleaning

### Training


## Details
### identification.py
package for identification of bacteria spectra
### parse_data.py
cleaning and transforming fasta file into csv file
### plotting.py
plotting fuctions in the paper
### training.py
training

## Publications
Ding Cheng, Liang Qiao and Peter Horvatovich, "Library-free MALDI-TOF MS Bacteria Typing", Submitted.

## License

BacteriaMS is distributed under a BSD license. See the LICENSE file for details.

## Contacts

Please report any problems directly to the github issue tracker. Also, you can send feedback to liang_qiao@fudan.edu.cn

