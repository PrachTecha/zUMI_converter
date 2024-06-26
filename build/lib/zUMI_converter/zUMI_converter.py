#!/usr/bin/env python
import os
import sys
import argparse
import pickle
from anndata import AnnData
from mudata import MuData
from scipy.sparse import csr_matrix
from pathlib import Path
# prevent rpy2 from using the R installed in the system
os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"


import rpy2.robjects as robjects
robjects.r['options'](warn=-1)
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import STAP
from rpy2.robjects.packages import importr
pandas2ri.activate()
importr('Matrix')

def zumis_output2mudata(input_path:str, output_path:str=None):
    """Convert RDS of zUMIs output to MuData format.

    Args:
        input_path (str): Path to RDS file of zUMIs output
        output_path (str, optional): Path to save the MuData object. Defaults to None.

    Returns:
        MuData: MuData object containing zUMIs output.
    """    
    mfunc = 'to_df <- function(dobj){return(as.data.frame(as.matrix(dobj)))}'
    rsparse2pandas = STAP(mfunc, "to_df")

    readRDS = robjects.r['readRDS']
    zumis_data = readRDS(input_path)
    zd = {k:v for k,v in zumis_data.items() if k in ['umicount', 'readcount']}
    for key in zd.keys():
        zd[key] = dict(zip(zd[key].names, list(zd[key])))
        for sub_key in ['exon', 'inex', 'intron']:
            zd[key][sub_key] = dict(zip(zd[key][sub_key].names, list(zd[key][sub_key])))
            zd[key][sub_key]['all'] = rsparse2pandas.to_df(zd[key][sub_key]['all'])
            zd[key][sub_key]['all'] = robjects.conversion.get_conversion().rpy2py(zd[key][sub_key]['all']).T
    
    adatas = {}
    for key in ['exon', 'inex', 'intron']:
        adata = AnnData(zd['umicount'][key]['all'], layers={'read': csr_matrix(zd['readcount'][key]['all'].values),
                                                             'umi': csr_matrix(zd['umicount'][key]['all'].values)})
        adata.X = csr_matrix(adata.X)
        adatas[key] = adata
    mdata = MuData(adatas)
    
    if output_path is not None:
        output_path = Path(output_path)
        if not output_path.parent.exists():
            print(f'Directory does not exists.\nCreating directory: {output_path.parent}')
            output_path.parent.mkdir(parents=True, exist_ok=True)
        mdata.write_h5mu(output_path, compression='gzip')
            
    return mdata

def argument_parser():
    parser = argparse.ArgumentParser(description="Convert RDS of zUMIs output to MuData")

    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument("-i", "--input-path", help="Path to RDS file of zUMIs output", type=str,
                          default=None, required=True)
    
    # Optional arguments
    parser.add_argument("-o", "--output-path", type=str, help="Path for MuData object to be stored",
                          default='./zUMI_counts.h5mu')

    
    arguments = parser.parse_args()

    return arguments.__dict__

def main():
    arguments = argument_parser()
    zumis_output2mudata(**arguments)

if __name__ == "__main__":
    main()