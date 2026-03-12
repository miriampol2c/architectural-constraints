import torch
import logging
import argparse
import numpy as np
from model import ProteinMPNN
from dataset import build_training_clusters, loader_pdb, get_pdbs, featurize, PDB_dataset, StructureDataset, StructureLoader
from utils import loss_nll, loss_smoothed

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
argparser.add_argument("--top_k", type=int, default=48, help="number of nearest neighbors to select")
argparser.add_argument("--data_path", type=str, default="pdb_2021aug02", help="path to the extracted dataset")
argparser.add_argument("--list_path", type=str, default="../original_dataset_filtered_by_all_enzymes.csv", help="path to csv file with cluster IDs")
argparser.add_argument("--model_path", type=str, default="model_weights/protein_mpnn_test_best.pt", help="path to pth file with model weights")

args = argparser.parse_args()
logging.basicConfig(filename='loggs/testing.log', level=logging.DEBUG)

data_path = args.data_path
params = {
    "LIST"    : args.list_path, 
    "VAL"     : f"{data_path}/valid_clusters.txt",
    "TEST"    : f"{data_path}/test_clusters.txt",
    "DIR"     : f"{data_path}",
    "DATCUT"  : "2030-Jan-01",
    "RESCUT"  : 3.5, #resolution cutoff for PDBs
    "HOMO"    : 0.70 #min seq.id. to detect homo chains
}


train, valid, test  = build_training_clusters(params, False)

test_dataset = PDB_dataset(list(test.keys()), loader_pdb, test, params)
logging.info(f"Test dataset clusters: {len(test_dataset)}")
test_pdb_list = get_pdbs(test_dataset)
test_structure_dataset = StructureDataset(test_pdb_list, max_length=1500)
test_loader = StructureLoader(test_structure_dataset)

top_k=args.top_k
temp = 0.1
node_hidden_dim = 128
device = torch.device("cuda")
checkpoint_path = args.model_path
checkpoint = torch.load(checkpoint_path, map_location=device)
model = ProteinMPNN(num_letters=21,
                        node_features=node_hidden_dim,
                        edge_features=node_hidden_dim,
                        hidden_dim=node_hidden_dim,
                        k_neighbors=top_k,
                        augment_eps=0.00)

model_parameters = filter(lambda p: p.requires_grad, model.parameters())
params = sum([np.prod(p.size()) for p in model_parameters])
logging.info(f"Model: {params} trainable parameters")
model.to(device)
model.load_state_dict(checkpoint['model_state_dict'])
model.eval()

with torch.no_grad():
    batch_lengths = []
    seq_rates = []
    test_sum, test_weights = 0.0, 0.0
    total_seq_req = 0.0
    ix = 0
    for batch in test_loader:
        X, S, mask, lengths, chain_M, residue_idx, mask_self, chain_encoding_all = featurize(batch, device)
        
        sample_dict = model.sample(X.to(device), S.to(device), chain_M.to(device), chain_encoding_all.to(device), residue_idx.to(device), mask=mask.to(device), temperature=temp)
        log_probs = model(X.to(device), S.to(device), mask.to(device), chain_M.to(device), residue_idx.to(device), chain_encoding_all.to(device), torch.randn(chain_M.shape, device=X.device))
        S_sample = sample_dict["S"]
            
        mask_for_loss = mask.to(device)*chain_M.to(device)
        loss, loss_av, true_false = loss_nll(S.to(device), log_probs, mask_for_loss.to(device))
        test_sum += torch.sum(loss * mask_for_loss.to(device)).cpu().data.numpy()
        test_weights += torch.sum(mask_for_loss).cpu().data.numpy()
        seq_recovery_rate = torch.sum(torch.sum(torch.nn.functional.one_hot(S.to(device), 21)*torch.nn.functional.one_hot(S_sample, 21),axis=-1)*mask_for_loss)/torch.sum(mask_for_loss)
        seq_rec_print = np.format_float_positional(np.float32(seq_recovery_rate.detach().cpu().numpy()), unique=False, precision=4)
        logging.info(f"Batch [{ix}] Sequence Recovery Rate: {seq_rec_print}")
        ix+=1
        batch_lengths.append(X.shape[0])
        seq_rates.append(np.float32(seq_recovery_rate.detach().cpu().numpy()))
        total_seq_req += torch.sum(torch.sum(torch.nn.functional.one_hot(S.to(device), 21)*torch.nn.functional.one_hot(S_sample, 21),axis=-1)*mask_for_loss)
        
    avg_seq_rate = np.average(np.array(seq_rates), weights=np.array(batch_lengths))
    test_loss = test_sum / test_weights
    perplexity = np.exp(test_loss)
    avg_seq_weighted_rate = total_seq_req / test_weights
    logging.info(f"Average Sequence Recovery Rate: {avg_seq_rate}")
    logging.info(f"Average Weighted Sequence Recovery Rate: {avg_seq_weighted_rate}")
    logging.info(f"Average perplexity: {perplexity}")
