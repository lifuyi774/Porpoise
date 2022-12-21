# Porpoise: computational analysis and prediction of RNA pseudouridine sites by a stacked machine learning framework.
## introduction
This study proposes a novel bioinformatics approach, termed Porpoise, for accurate identifying RNA pseudouridine sites. Porpoise is developed by comprehensively evaluated 18 popular feature encoding schemes and four types of features, including binary features, pseudo k-tuple composi-tion (PseKNC), nucleotide chemical property (NCP) and position-specific trinucleotide propensity based on single-strand (PSTNPss) are selected and feed into the stacked framework to construct-ing the final meta-learning model.Cross-validation tests on the benchmark dataset and independent tests demonstrate that Porpoise achieved superior predictive performance than state-of-the-art approaches.

This is a repository of codes of Porpoise, an RNA pseudouridine sites predictor. You can use this program and know more about it through our [website]
(http://web.unimelb-bioinfortools.cloud.edu.au/Porpoise/).

## Environment

The environment on our computer is as follows:
*python 3.5~3.6
* mlxtend 0.13.0
* Pandas 0.24.2
* NumPy 1.18.5
* scikit-learn 0.19.0
*xgboost 0.90

## Data

-RNA_training.txt: The benchmark dataset S(1) for H. sapiens. It is formed by 495 Ψ-site-containing sequences and 495 false Ψ-site-containing sequences. Each of these samples is 21-bp long with the uridine located at the center. None of the sequences included here has ≥60% pairwise sequence identity to any other in a same subset. 
-RNA_training1.txt: The benchmark dataset S(2) for S. cerevisiae. It is formed by 314 Ψ-site-containing sequences and 314 false ψ-site-containing sequences. Each of these samples is 31-bp long with the uridine located at the center. None of the sequences included here has ≥60% pairwise sequence identity to any other in a same subset.
-RNA_training2.txt: The benchmark dataset S(3) for M. musculus. It is formed by 472 Ψ-site-containing sequences and 472 false ψ-site-containing sequences. Each of these samples is 21-bp long with the uridine located at the center. None of the sequences included here has ≥60% pairwise sequence identity to any other in a same subset.
- RNA_test.txt: The independent dataset S(4) for H. sapiens. It is formed by 100 ψ-site-containing sequences and 100 false ψ-site-containing sequences. Each sample is 21-bp long with uridine located at the center. None of the samples included here occurs in S(1) of the Supporting Information S1.
- RNA_test1.txt: The independent dataset S(5) for S. cerevisiae. It is formed by 100 ψ-site-containing sequences and 100 false ψ-site-containing sequences. Each sample is 31-bp long with the uridine located at the center. None of the samples included here occurs in S(2) of the Supporting Information S2.

 These  benchmark datasets of pseudouridylation sites were taken from the additional materials of Chen et al.(Chen, W., Tang, H., Ye, J., Lin, H., and Chou, K.-C. (2016). iRNA-PseU: Identifying
RNA pseudouridine sites. Mol. Ther. Nucleic Acids 5, e332.)

## Usage

The parameters of the model can be passed as arguments to the script. 
It is possible to see which parameters can be passed by running (from inside Porpoise/):

```
cd Porpoise/iuse
H990_Model.py/M944_Model.py/S628_Model.py/pre-processing.py/MYSVG.py  -h
```

Use the following  command to predict (from inside Porpoise/iuse/):

```
python H990_Model.py --i='path/to/your/fasta file' --o='path/to/your/result.txt'
python M944_Model.py --i='path/to/your/fasta file' --o='path/to/your/result.txt'
python S628_Model.py --i='path/to/your/fasta file' --o='path/to/your/result.txt'
```
Check the length of the input sequence (The minimum sequence length of H. sapiens and M. musculus is 21 and  the minimum sequence length of S. cerevisiae is 31)


Use the following command to create your stack model (from inside Porpoise/iuse/):
```
./Auto-pipeline --i='path/to/your/fasta file' --c config.txt --o='path/to/your/result.txt'
```
Use the following command to create the SVG (from inside Porpoise/iuse/):
```
python MYSVG.py --i='path/to/your/fasta file' --t threshold --m ='path/to/your/model/result.txt' --s species --o ='path/to/your/model/result.svg'
```