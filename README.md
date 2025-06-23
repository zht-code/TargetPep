# TargetPep: Multi-attribute Prompt-based Large Language Model Fine-Tuning for Target-Aware Peptide Design
<p align="center">
    <img src="temp/img.png" width="800" class="center" alt="TargetPep Workflow"/>
    <br/>
</p>
## Installation

#### Create the conda environment and activate it.
```
conda create -n TargetPep python==3.10.14
conda activate TargetPep
```
#### Install basic packages
```
# install requirements
pip install -r requirements.txt

pip install easydict
pip install biopython
# mmseq
conda install bioconda::mmseqs2

# Alternative: obabel and RDkit
conda install -c openbabel openbabel
conda install conda-forge::rdkit

# Alternative for visualization: py3dmol
conda install conda-forge::py3dmol

#Install the sequence similarity comparison tool EMBOSS
conda install -c bioconda emboss
```

### Packages for training and generating.

#### Install pytorch 2.4.1 with the cuda version that is compatible with your device.
```
# torch-scatter
pip install torch-scatter -f https://data.pyg.org/whl/torch-2.4.1+cu118.html  

```

## Dataset 
We provide the dataset of `PPBench2024` through [google drive](https://drive.google.com/drive/folders/1Fcxh45AdlFfpimFZEyK_lKOBq2HmEcH2?usp=drive_link), together with `PPDBench'.

Please download `PPBench2024.zip` `PPDBench.zip` and unzip it, leading to the data file directory as 
```
- dataset
    - Train
        - 1a0m_A
            peptide.pdb
            recepotor.pdb
    - Test
        - 1cjr
            peptide.pdb
            recepotor.pdb
        - 1cka
            peptide.pdb
            recepotor.pdb
        ...      
```

## Training and Generating
### Training from scratch
Run the following command for TargetPep training:

```
python train.py
```

#### After training, you can choose an epoch for generating the peptides through:

For batch generation, you should first download the pre-trained `best_model_epoch_25_loss_0.4112.pth` model from [google drive](https://drive.google.com/drive/folders/1Fcxh45AdlFfpimFZEyK_lKOBq2HmEcH2?usp=drive_link), then save it as `. /checkpoint/best_model_epoch_25_loss_0.4112.pth` and finally run the following command for batch generation:
For ESM3 and Progen2 weights, you can download the `log.zip` configuration file from [google drive](https://drive.google.com/file/d/17Q8rNNs97jdl1F_MgEPY7JDQhZpiAknK/view?usp=sharing).
```
python ./utils/evaluate/batch_generate.py

python ./utils/evaluate/batch_generate_top3.py

python ./utils/evaluate/batch_generate_top5.py
```
If you want to evaluate peptides directly, they are available from our Google Drive as `peptides1-5.zip`, containing 5 samples per protein structure for more robust evaluation.


## Packages and Scripts for Evaluation

### Packages for docking and other evaluation.

#### Vina: For Vina Docking, install the packages through:
```
 conda install conda-forge::vina
 pip install meeko
 pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3.git@aee55d50d5bdcfdbcd80220499df8cde2a8f4b2a
 pip install pdb2pqr
```
`./utils/datasets/vina_score_computing.py` gives an example of our Python interface for vinadock.

#### HDock: For HDock, firstly, libfftw3 is needed for hdock with `apt-get install -y libfftw3-3`. Besides, the HDock software can be downloaded through: http://huanglab.phys.hust.edu.cn/software/hdocklite/. After downloading it, install or unzip it to the `./bin` directory, leading to the file structure as 
```
- bin
    - hdock
        1CGl_l_b.pdb
        1CGl_r_b.pdb
        createpl
        hdock.out
```
`./utils/evaluate/hdock.py`  gives an example of our python interface for hdock.

`./utils/evaluate/hdock_PPflow.py`  gives an example of PPflow python interface for hdock.

`./utils/evaluate/hdock_RFdiffusion.py`  gives an example of RFdiffusion python interface for hdock.

#### TMscore: The available TMscore evaluation software is provided in `./data/data/downloads`, as 
```
- data
	-data
		-downloads
			TMalign
			TMalign.cpp
```

### Evaluation scripts
The evaluation scripts are given in `./utils/evaluate` directory, you can run the following for the evaluation in the main experiments:

```
# Evaluting the physical chemical properties
python ./utils/evaluate/physical_chemical_properties.py

# Evaluting the Diversity
./utils/evaluate/diversity.ipynb

# Evaluting the Solubility
./utils/evaluate/solubility.py

# Evaluting the sc-TM
./utils/evaluate/scTM.ipynb

# Evaluting the the-state-of-the-art
./utils/evaluate/ablationquvina.py
./utils/evaluate/ablation_solubility.py
```

