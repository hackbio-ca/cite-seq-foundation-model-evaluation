from mudata import MuData

import numpy as np 
import anndata as ann 
import mudata as md 
from sklearn.model_selection import train_test_split
import torch
import pandas as pd


def scenario_1_split(mudata_object, train_size = 0.6, val_size = 0.1, test_size = 0.3, seed = 42):
    """
    Split the MuonData object into training, validation, and test sets. In this scenario, the split is 
    based on the rows (cells) of the MuonData object.
    
    Parameters:
    mudata_object (MuonData): MuonData object to split
    train_size (float): Proportion of the dataset to include in the training set
    val_size (float): Proportion of the dataset to include in the validation set
    test_size (float): Proportion of the dataset to include in the test set
    
    Returns:
    train_data (MuonData): MuonData object with annotations for train, val, test splits
    """
    # Get the indices of the mdata by enumerating over shape
    mudata_indices = [i for i in range(mudata_object.shape[0])]
    
    # Set the seed for reproducibility
    np.random.seed(seed)
    
    # Split the indices based on the specified test/train/validation sizes
    train_indices, test_indices = train_test_split(mudata_indices, train_size=train_size)
    val_indices, test_indices = train_test_split(test_indices, train_size=(val_size/(val_size + test_size)))
    
    # Add the indices to the obs of the ADT and SCT aspects of the MuonData object
    mudata_object["ADT"].obs["Split"] = "None"
    mudata_object["ADT"].obs["Split"][train_indices] = "Train"
    mudata_object["ADT"].obs["Split"][val_indices] = "Validation"
    mudata_object["ADT"].obs["Split"][test_indices] = "Test"
    
    mudata_object["SCT"].obs["Split"] = "None"
    mudata_object["SCT"].obs["Split"][train_indices] = "Train"
    mudata_object["SCT"].obs["Split"][val_indices] = "Validation"
    mudata_object["SCT"].obs["Split"][test_indices] = "Test"
    
    # Return the annotated 
    return mudata_object


def scenario_2_split(mudata_object, train_size = 0.7, val_size = 0.15, test_size = 0.15, seed = 42):
    """
    Split the MuonData object into training, validation, and test sets. In this case, the split is 
    amongst the vars of the ADT aspect of the MuonData object.
    
    Parameters:
    mudata_object (MuonData): MuonData object to split
    train_size (float): Proportion of the dataset to include in the training set
    val_size (float): Proportion of the dataset to include in the validation set
    test_size (float): Proportion of the dataset to include in the test set
    
    Returns:
    mudata_object (MuonData): MuonData object with the ADT vars annotated with the split
    """
    # Get the indices of the ADT vars by enumerating over shape
    adt_var_indices = [i for i in range(mudata_object["ADT"].var.shape[0])]
    
    # Set the seed for reproducibility
    np.random.seed(seed)
    
    # Split the indices based on the specified test/train/validation sizes
    train_indices, test_indices = train_test_split(adt_var_indices, train_size=train_size)
    val_indices, test_indices = train_test_split(test_indices, train_size=(val_size/(val_size + test_size)))
    
    
    # Subset the mudata objects based on indices and return train/test/validation
    mudata_object["ADT"].var["Split"] = "None"
    mudata_object["ADT"].var["Split"][train_indices] = "Train"
    mudata_object["ADT"].var["Split"][val_indices] = "Validation"
    mudata_object["ADT"].var["Split"][test_indices] = "Test"
    
    # Return the annodated split data
    return mudata_object

def adt_split_returns(mudata_train, scenario=1, return_type="mudata"):
    """
    Function to take mudata object returns from Scenario 1 and 2 splitters, and return either anndata
    slices of ADT and SCT, numpy arrays of ADT and SCT, or torch tensors of ADT and SCT.
    
    Args:
        scenario (int, optional): Type of scenario. Defaults to 1. 1 or 2. 
        mudata_train (_type_, optional): Train mudata object. Defaults to None
        return_type (str, optional): Either "AnnData", "Numpy", or "Torch. Defaults to "AnnData"

    Returns:
        Dependant on return_type, returns either anndata, numpy, or torch tensors of ADT and SCT data.
    """
    # Scenario 1 conditions 
    if scenario == 1:
        if return_type == "AnnData":
            train_adt_adata = mudata_train["ADT"][mudata_train["ADT"].obs["Split"] == "Train"].copy()
            train_sct_adata = mudata_train["SCT"][mudata_train["SCT"].obs["Split"] == "Train"].copy()
            
            val_adt_adata = mudata_train["ADT"][mudata_train["ADT"].obs["Split"] == "Validation"].copy()
            val_sct_adata = mudata_train["SCT"][mudata_train["SCT"].obs["Split"] == "Validation"].copy()
            
            test_adt_adata = mudata_train["ADT"][mudata_train["ADT"].obs["Split"] == "Test"].copy()
            test_sct_adata = mudata_train["SCT"][mudata_train["SCT"].obs["Split"] == "Test"].copy()
            
            return train_adt_adata, train_sct_adata, val_adt_adata, val_sct_adata, test_adt_adata, test_sct_adata

        elif return_type == "mudata":
            train_adt_adata = mudata_train["ADT"][mudata_train["ADT"].obs["Split"] == "Train"].copy()
            train_sct_adata = mudata_train["SCT"][mudata_train["SCT"].obs["Split"] == "Train"].copy()
            
            val_adt_adata = mudata_train["ADT"][mudata_train["ADT"].obs["Split"] == "Validation"].copy()
            val_sct_adata = mudata_train["SCT"][mudata_train["SCT"].obs["Split"] == "Validation"].copy()
            
            test_adt_adata = mudata_train["ADT"][mudata_train["ADT"].obs["Split"] == "Test"].copy()
            test_sct_adata = mudata_train["SCT"][mudata_train["SCT"].obs["Split"] == "Test"].copy()
            
            train = MuData({"ADT": train_adt_adata, "SCT": train_sct_adata})
            val = MuData({"ADT": val_adt_adata, "SCT": val_sct_adata})
            test = MuData({"ADT": test_adt_adata, "SCT": test_sct_adata})
            return train, val, test
        
        elif return_type == "Numpy":
            train_adt_adata = mudata_train["ADT"][mudata_train["ADT"].obs["Split"] == "Train"].X
            train_sct_adata = mudata_train["SCT"][mudata_train["SCT"].obs["Split"] == "Train"].X.todense()
            
            val_adt_adata = mudata_train["ADT"][mudata_train["ADT"].obs["Split"] == "Validation"].X
            val_sct_adata = mudata_train["SCT"][mudata_train["SCT"].obs["Split"] == "Validation"].X.todense()
            
            test_adt_adata = mudata_train["ADT"][mudata_train["ADT"].obs["Split"] == "Test"].X
            test_sct_adata = mudata_train["SCT"][mudata_train["SCT"].obs["Split"] == "Test"].X.todense()
            
            return train_adt_adata, train_sct_adata, val_adt_adata, val_sct_adata, test_adt_adata, test_sct_adata
            
        elif return_type == "Torch":
            train_adt_adata = torch.tensor(mudata_train["ADT"][mudata_train["ADT"].obs["Split"] == "Train"].X)
            train_sct_adata = torch.tensor(mudata_train["SCT"][mudata_train["SCT"].obs["Split"] == "Train"].X.todense())
            
            val_adt_adata = torch.tensor(mudata_train["ADT"][mudata_train["ADT"].obs["Split"] == "Validation"].X)
            val_sct_adata = torch.tensor(mudata_train["SCT"][mudata_train["SCT"].obs["Split"] == "Validation"].X.todense())
            
            test_adt_adata = torch.tensor(mudata_train["ADT"][mudata_train["ADT"].obs["Split"] == "Test"].X)
            test_sct_adata = torch.tensor(mudata_train["SCT"][mudata_train["SCT"].obs["Split"] == "Test"].X.todense())
            
            return train_adt_adata, train_sct_adata, val_adt_adata, val_sct_adata, test_adt_adata, test_sct_adata
        else:
            raise ValueError("Return type must be 'AnnData', 'Numpy', or 'Torch'")
        
    # Scenario 2 conditions 
    elif scenario == 2:
        if return_type == "AnnData":
            train_sct_adata = mudata_train["SCT"].copy()
            val_sct_adata = mudata_train["SCT"].copy()
            test_sct_adata = mudata_train["SCT"].copy()
            
            train_adt_adata = mudata_train["ADT"][:, mudata_train["ADT"].var["Split"] == "Train"].copy()
            val_adt_adata = mudata_train["ADT"][:, mudata_train["ADT"].var["Split"] == "Validation"].copy()
            test_adt_adata = mudata_train["ADT"][:, mudata_train["ADT"].var["Split"] == "Test"].copy()
            
            return train_adt_adata, train_sct_adata, val_adt_adata, val_sct_adata, test_adt_adata, test_sct_adata

        elif return_type == "mudata":
            train_sct_adata = mudata_train["SCT"].copy()
            val_sct_adata = mudata_train["SCT"].copy()
            test_sct_adata = mudata_train["SCT"].copy()
            
            train_adt_adata = mudata_train["ADT"][:, mudata_train["ADT"].var["Split"] == "Train"].copy()
            val_adt_adata = mudata_train["ADT"][:, mudata_train["ADT"].var["Split"] == "Validation"].copy()
            test_adt_adata = mudata_train["ADT"][:, mudata_train["ADT"].var["Split"] == "Test"].copy()
            
            train = MuData({"ADT": train_adt_adata, "SCT": train_sct_adata})
            val = MuData({"ADT": val_adt_adata, "SCT": val_sct_adata})
            test = MuData({"ADT": test_adt_adata, "SCT": test_sct_adata})
            return train, val, test

        elif return_type == "Numpy":
            train_sct_adata = mudata_train["SCT"].X.todense()
            val_sct_adata = mudata_train["SCT"].X.todense()
            test_sct_adata = mudata_train["SCT"].X.todense()
            
            train_adt_adata = mudata_train["ADT"][:, mudata_train["ADT"].var["Split"] == "Train"].X
            val_adt_adata = mudata_train["ADT"][:, mudata_train["ADT"].var["Split"] == "Validation"].X
            test_adt_adata = mudata_train["ADT"][:, mudata_train["ADT"].var["Split"] == "Test"].X
            
            return train_adt_adata, train_sct_adata, val_adt_adata, val_sct_adata, test_adt_adata, test_sct_adata
            
        elif return_type == "Torch":
            train_sct_adata = torch.tensor(mudata_train["SCT"].X.todense())
            val_sct_adata = torch.tensor(mudata_train["SCT"].X.todense())
            test_sct_adata = torch.tensor(mudata_train["SCT"].X.todense())
            
            train_adt_adata = torch.tensor(mudata_train["ADT"][:, mudata_train["ADT"].var["Split"] == "Train"].X)
            val_adt_adata = torch.tensor(mudata_train["ADT"][:, mudata_train["ADT"].var["Split"] == "Validation"].X)
            test_adt_adata = torch.tensor(mudata_train["ADT"][:, mudata_train["ADT"].var["Split"] == "Test"].X)
            
            return train_adt_adata, train_sct_adata, val_adt_adata, val_sct_adata, test_adt_adata, test_sct_adata
            
        else:
            raise ValueError("Return type must be 'AnnData', 'Numpy', or 'Torch'")
        
        
def load_appris(unique_transcripts=True):
    # generate doc string
    """
    Load the appris data
    :param unique_transcripts: whether to load only unique transcripts
    :return: the appris data
    """
    # ## load human appris
    dir = '/h/phil/Documents/01_projects/contrastive_rna_representation/'

    app_h = pd.read_csv(f'{dir}/data/appris_data_human.principal.txt', sep='\t')
    print(app_h['Gene ID'].duplicated().sum())
    app_h['numeric_value'] = app_h['APPRIS Annotation'].str.split(':').str[1]
    app_h['key_value'] = app_h['APPRIS Annotation'].str.split(':').str[0]
    app_h = app_h.sort_values(
        ['Gene ID', 'key_value', 'numeric_value', "Transcript ID"],
        ascending=[True, False, True, True],
    )
    if unique_transcripts:
        app_h = app_h[~app_h.duplicated('Gene ID')]
        app_h = app_h[~app_h.duplicated('Gene name')]
    return app_h


def data_loader_wrapper(mdata, split_type='genes', seed=42):
    assert split_type in ['genes', 'cells', 'genes_cells'], "Split type must be either 'genes' or 'cells'"

    print('Before Filtering')
    print(mdata['ADT'].X.shape, mdata['ADT'].var.shape)
    mdata2 = mdata['ADT'][:, mdata["ADT"].var['gene_name'].notnull().values]
    blacklist_genes = [
        'SIGLEC8', 'SELE', 'CDH17','THY1', 
        'CD177', 'KDR', 'CEACAM8', 'VTCN1', 'TEK'
    ]
    mdata2 = mdata2[:, ~mdata2.var['gene_name'].isin(blacklist_genes)]

    # drop rows with duplicate gene names drop both
    mdata2 = mdata2[:, ~mdata2.var.duplicated('gene_name', keep=False)]
    
    # drop rows with esm_0 null
    mdata2 = mdata2[:, mdata2.var['esm_0'].notnull()]
    mdata2 = mdata2[:, mdata2.var['orthrus_0'].notnull()]
    
    mdata = MuData({"ADT": mdata2, "SCT": mdata['SCT']})
    print('After Filtering')
    print(mdata['ADT'].X.shape, mdata['ADT'].var.shape)
    
    if split_type == 'genes':
        mdata = scenario_2_split(mdata, seed=seed)
        scenario = 2
        
    elif split_type == 'cells':
        mdata = scenario_1_split(mdata, seed=seed)
        scenario = 1
    
    elif split_type == 'genes_cells':
        raise NotImplementedError("Split type 'genes_cells' not implemented yet")
    
    train, val, test = adt_split_returns(scenario=scenario, mudata_train = mdata, return_type = "mudata")
    return train, val, test
    
    