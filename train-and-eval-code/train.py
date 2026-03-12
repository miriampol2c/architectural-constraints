import torch
import logging
import argparse
import numpy as np
from model import ProteinMPNN
from dataset import build_training_clusters, loader_pdb, get_pdbs, featurize, PDB_dataset, StructureDataset, StructureLoader
from utils import loss_nll, loss_smoothed, get_std_opt

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
argparser.add_argument("--data_path", type=str, default="pdb_2021aug02", help="path to the extracted dataset")
argparser.add_argument("--list_path", type=str, default="../original_dataset_filtered_by_all_enzymes.csv", help="path to csv file with cluster IDs")
argparser.add_argument("--epochs", type=int, default=150, help="number of epochs to train for")
argparser.add_argument("--top_k", type=int, default=48, help="number of nearest neighbors to select")
argparser.add_argument("--reload_data_every_n_epochs", type=int, default=2, help="reload training data every n epochs")
argparser.add_argument("--backbone_noise", type=float, default=0.02, help="amount of noise added to backbone during training") 

args = argparser.parse_args()
logging.basicConfig(filename='loggs/training.log', level=logging.DEBUG)

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

train_dataset = PDB_dataset(list(train.keys()), loader_pdb, train, params)
logging.info(f"Train dataset clusters: {len(train_dataset)}")
train_pdb_list = get_pdbs(train_dataset)
train_structure_dataset = StructureDataset(train_pdb_list, max_length=1500)
train_loader = StructureLoader(train_structure_dataset)

valid_dataset = PDB_dataset(list(valid.keys()), loader_pdb, valid, params)
logging.info(f"Valid dataset clusters: {len(valid_dataset)}")
valid_pdb_list = get_pdbs(valid_dataset)
valid_structure_dataset = StructureDataset(valid_pdb_list, max_length=1500)
valid_loader = StructureLoader(valid_structure_dataset)

top_k=args.top_k
node_hidden_dim = 128
DEVICE = torch.device("cuda")
epochs = args.epochs
model = ProteinMPNN(num_letters=21,
                        node_features=node_hidden_dim,
                        edge_features=node_hidden_dim,
                        hidden_dim=node_hidden_dim,
                        k_neighbors=top_k,
                        augment_eps=args.backbone_noise)

model.to(DEVICE)
total_step = 0
optimizer = get_std_opt(model.parameters(), node_hidden_dim, total_step)
best_loss = 100000000
reload_c = 0

for epoch in range(epochs):
    if epoch % args.reload_data_every_n_epochs == 0:
        if reload_c != 0:
            train_dataset = PDB_dataset(list(train.keys()), loader_pdb, train, params)
            logging.info(f"Train dataset clusters: {len(train_dataset)}")
            train_pdb_list = get_pdbs(train_dataset)
            train_structure_dataset = StructureDataset(train_pdb_list, max_length=1500)
            train_loader = StructureLoader(train_structure_dataset)

            valid_dataset = PDB_dataset(list(valid.keys()), loader_pdb, valid, params)
            logging.info(f"Valid dataset clusters: {len(valid_dataset)}")
            valid_pdb_list = get_pdbs(valid_dataset)
            valid_structure_dataset = StructureDataset(valid_pdb_list, max_length=1500)
            valid_loader = StructureLoader(valid_structure_dataset)
        reload_c += 1
        
    epoch_loss = 0.0
    train_acc = 0.0
    train_sum, train_weights = 0., 0.
    num_samples = 0
    model.train()
    for batch in train_loader:
        X, S, mask, lengths, chain_M, residue_idx, mask_self, chain_encoding_all = featurize(batch, DEVICE)
        optimizer.zero_grad()
        log_probs = model(X.to(DEVICE), S.to(DEVICE), mask.to(DEVICE), chain_M.to(DEVICE), residue_idx.to(DEVICE), chain_encoding_all.to(DEVICE), torch.randn(chain_M.shape, device=X.device))
        mask_for_loss = mask*chain_M
        
        _, loss_av_smoothed = loss_smoothed(S.to(DEVICE), log_probs, mask_for_loss.to(DEVICE))
        loss_av_smoothed.backward()

        #if args.gradient_norm > 0.0:
        #    total_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), args.gradient_norm)

        optimizer.step()
            
        loss, loss_av, true_false = loss_nll(S.to(DEVICE), log_probs, mask_for_loss.to(DEVICE))
            
        train_sum += torch.sum(loss * mask_for_loss.to(DEVICE)).cpu().data.numpy()
        train_acc += torch.sum(true_false * mask_for_loss.to(DEVICE)).cpu().data.numpy()
        train_weights += torch.sum(mask_for_loss).cpu().data.numpy()

        total_step += 1
        
    train_loss = train_sum / train_weights
    train_accuracy = train_acc / train_weights
    train_perplexity = np.exp(train_loss)
    
    logging.info(f"[{epoch}] Train Loss: {train_loss}")
    logging.info(f"[{epoch}] Train Accuracy: {train_accuracy}, Train Perplexity: {train_perplexity}")
    
    model.eval()
    with torch.no_grad():
        validation_acc = 0.0
        validation_sum, validation_weights = 0., 0.
        num_val_samples = 0
        for batch in valid_loader:
            X, S, mask, lengths, chain_M, residue_idx, mask_self, chain_encoding_all = featurize(batch, DEVICE)
            log_probs = model(X.to(DEVICE), S.to(DEVICE), mask.to(DEVICE), chain_M.to(DEVICE), residue_idx.to(DEVICE), chain_encoding_all.to(DEVICE), torch.randn(chain_M.shape, device=X.device))
            mask_for_loss = mask*chain_M
            loss, loss_av, true_false = loss_nll(S.to(DEVICE), log_probs, mask_for_loss.to(DEVICE))
                    
            validation_sum += torch.sum(loss * mask_for_loss.to(DEVICE)).cpu().data.numpy()
            validation_acc += torch.sum(true_false * mask_for_loss.to(DEVICE)).cpu().data.numpy()
            validation_weights += torch.sum(mask_for_loss).cpu().data.numpy()
                
    validation_loss = validation_sum / validation_weights
    validation_accuracy = validation_acc / validation_weights
    validation_perplexity = np.exp(validation_loss)
    
    logging.info(f"[{epoch}] Valid Loss: {validation_loss}")
    logging.info(f"[{epoch}] Valid Accuracy: {validation_accuracy}, Valid Perplexity: {validation_perplexity}")
    
    if validation_loss < best_loss:
        best_loss = validation_loss
        torch.save({'model_state_dict': model.state_dict()}, f'model_weights/protein_mpnn_test_best.pt')
        logging.info("Saved best performing model...")
    
    torch.save({'model_state_dict': model.state_dict()}, f'model_weights/protein_mpnn_test_last.pt')
    logging.info("*******************")
