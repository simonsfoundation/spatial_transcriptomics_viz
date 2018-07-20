# spatial_transcriptomics_viz - visualization aids for spatial transcriptomics

## Setup instructions

You will need a jupyter installation.  I recommend
installing anaconda to make everything easier.

```
mkdir ~/.jupyter
mkdir ~/.jupyter/custom
cp custom.css ~/.jupyter/custom
conda create -n ipykernel_py2 python=2 ipykernel
source activate ipykernel_py2
python -m ipykernel install --user
conda install PIL
pip install -r requirements.txt
python setup.py install
jupyter nbextension enable --py --sys-prefix widgetsnbextension
```
