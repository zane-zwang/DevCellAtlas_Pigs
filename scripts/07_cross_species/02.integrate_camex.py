# Purpose: integrate camex.
import os
from pathlib import Path
PUBLICATION_ROOT = Path(os.environ.get("DEVPIGATLAS_ROOT", ".")).resolve()

import time
import torch
import shutil
import warnings
import argparse
import importlib
import scanpy as sc

import pandas as pd
import numpy as np

from CAMEX.base import Dataset
from CAMEX.trainer import Trainer

from runpy import run_path
PARAMS = run_path(str(Path(__file__).with_name('01.camex_params.py')))['PARAMS']

t1 = time.time()

time_start = time.strftime("%Y-%m-%d-%H-%M-%S")
log_path = f'./log/{time_start}/'
for k, v in PARAMS.items():
    v['time_start'] = time_start
    v['log_path'] = log_path
print(log_path)

os.makedirs(log_path, exist_ok=True)
shutil.copy(str(PUBLICATION_ROOT / 'scripts/07_cross_species/01.camex_params.py'), log_path + 'params_current.py')
print(f'time: {time_start}')

print('start preprocess')
dataset = Dataset(**PARAMS['preprocess'])
adata_CAMEX = dataset.adata_whole
dgl_data = dataset.dgl_data

print('start train')
trainer = Trainer(adata_CAMEX, dgl_data, **PARAMS['train'])

trainer.integration()

adata_CAMEX.write_h5ad(log_path + 'adata_CAMEX.h5ad', compression='gzip')
t2 = time.time()
print(f'time usage: {round(t2-t1)} seconds')
