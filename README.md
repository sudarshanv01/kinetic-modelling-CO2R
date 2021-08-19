Unified Mechanistic Understanding of CO2 Reduction to CO on Transition Metal and Single Atom Catalysts
-------------------------

Contains all data required to reproduce figures in [this](10.26434/chemrxiv.14427986) manuscript.

## Instructions

Create a new virtual environment with Python, this can be done as follows

```
python3 -m venv <environment_name>
source <environment_name>/bin/activate
```

To get all the required packages that were used to make the figures, 

```
python -m pip install -r requirements.txt
```

In order to get the AiiDA nodes for the CatMAP calculations, also do 

```
python -m pip install -r optional-requirements.txt
```

Instructions on how to reproduce each of the main figures are in `kinetic_modelling`