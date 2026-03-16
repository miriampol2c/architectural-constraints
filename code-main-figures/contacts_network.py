import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import numpy as np
import os
from collections import Counter
import random
import matplotlib.patches as mpatches

def read_and_filter(file, residues, label):
    df = pd.read_csv(file, sep='\t')
    df['MapRes1ID'] = df['NumRes1_Ref'].astype(str) + df['AA1']
    df['MapRes2ID'] = df['NumRes2_Ref'].astype(str) + df['AA2']

    filt1 = df[df['NumRes1_Ref'].isin(residues)]
    filt2 = df[df['NumRes2_Ref'].isin(residues)]
    filtered = pd.concat([filt1, filt2]).drop_duplicates()

    conserved = filtered[(filtered['FreqConts'] >= 0.5) & (filtered['ICtotal'] >= 0.5)]
    conserved['PairKey'] = conserved.apply(lambda row: tuple(sorted((row['NumRes1_Ref'], row['NumRes2_Ref']))), axis=1)
    conserved.to_csv(f"filtered_contacts_{label}.tsv", sep='\t', index=False)
    
    return conserved

def plot_barplot(conserved_list, labels, output_prefix):
    combined = pd.concat(conserved_list, keys=labels, names=['Dataset'])
    summary = combined.groupby(['Dataset', 'FstConserved']).size().reset_index(name='Count')

    sns.set(style='white', context='paper')
    palette = {'MIN': 'limegreen', 'NEU': 'darkgrey', 'MAX': 'red'}

    plt.figure(figsize=(6.5, 4.5))  # Optimized for single-column layout
    ax = sns.barplot(data=summary, x='Dataset', y='Count', hue='FstConserved',
                     palette=palette, edgecolor='black')

    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.annotate(f'{int(height)}', (p.get_x() + p.get_width() / 2., height),
                        ha='center', va='bottom', fontsize=8, xytext=(0, 3),
                        textcoords='offset points')

    ax.set_title("", fontsize=12, weight='bold')
    ax.set_ylabel("N conserved contacts involving catalytic residues", fontsize=10, weight='bold')
    ax.set_xlabel("")
    ax.legend(title="FrustIC", loc='upper right', frameon=False)

    # Tick formatting
    plt.xticks(fontsize=10, rotation=0)
    plt.yticks(fontsize=10)
    sns.despine()  # Clean look

    plt.tight_layout(pad=0.5)
    plt.savefig(f"{output_prefix}_barplot.svg", dpi=600, bbox_inches='tight')
    #plt.savefig(f"{output_prefix}_barplot.png", dpi=600, bbox_inches='tight')
    plt.close()
    print(f"Barplot saved as {output_prefix}_barplot.svg and .png")

import matplotlib.patches as mpatches

def plot_barplot_with_shared(conserved_list, labels, shared_keys, output_prefix):
    combined = pd.concat(conserved_list, keys=labels, names=['Dataset']).reset_index(level='Dataset')
    combined['IsShared'] = combined['PairKey'].apply(lambda x: x in shared_keys)

    # Total bars
    total_counts = combined.groupby(['Dataset', 'FstConserved']).size().reset_index(name='Count')

    # Shared portion only (same categories)
    shared_counts = (
        combined[combined['IsShared']]
        .groupby(['Dataset', 'FstConserved'])
        .size()
        .reset_index(name='SharedCount')
    )

    merged = pd.merge(total_counts, shared_counts, on=['Dataset', 'FstConserved'], how='left').fillna(0)
    merged['SharedCount'] = merged['SharedCount'].astype(int)

    sns.set(style='white', context='paper')
    palette = {'MIN': 'limegreen', 'NEU': 'darkgrey', 'MAX': 'red'}

    fig, ax = plt.subplots(figsize=(11.5, 7.5))

    # Plot total bars
    base = sns.barplot(data=merged, x='Dataset', y='Count', hue='FstConserved',
                       palette=palette, edgecolor='black', ax=ax)

    # Overlay hatched bars for shared contacts
    for i, bar in enumerate(base.patches):
        # Extract corresponding shared value
        group = merged.iloc[i // 3]  # Each dataset*FrustIC group has 3 bars
        shared_val = group['SharedCount']
        if shared_val > 0:
            x = bar.get_x()
            width = bar.get_width()
            height = shared_val
            bottom = bar.get_y()

            # Overlay hatched rectangle
            ax.add_patch(plt.Rectangle(
                (x, bottom), width, height,
                hatch='//', fill=False, edgecolor='black', linewidth=0, zorder=7
            ))

        # Annotate total bar height
        height = bar.get_height()
        if height > 0:
            ax.annotate(f'{int(height)}', (bar.get_x() + width / 2., height),
                        ha='center', va='bottom', xytext=(0, 3),
                        textcoords='offset points', fontsize=16)

    ax.set_ylabel("N conserved contacts involving catalytic residues", fontsize=18, weight='bold')
    ax.set_xlabel("")
    
    # Get current legend handles and labels from seaborn barplot
    handles, labels_ = ax.get_legend_handles_labels()

    # Create custom hatch legend handle
    hatch_handle = mpatches.Patch(facecolor='white', edgecolor='black', hatch='//', label='Shared contacts')

    # Append the hatch legend handle to the existing legend
    handles.append(hatch_handle)
    labels_.append('Shared contacts')

    # Draw updated legend
    ax.legend(handles=handles, labels=labels_, title="FrustIC", loc='upper right', frameon=False, fontsize=16, title_fontsize=16)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    sns.despine()
    plt.tight_layout(pad=0.5)

    plt.savefig(f"{output_prefix}_barplot_shared_overlay.svg", dpi=600, bbox_inches='tight')
    plt.close()
    print(f"Overlay barplot with shared markers saved as {output_prefix}_barplot_shared_overlay.svg")



def assign_colors(residues_of_interest):
    accent = list(plt.get_cmap('Accent').colors)
    manual_map = {
        70: accent[1],   # Purple
        73: accent[3],   # Yellow
        130: accent[0],  # Green
        166: accent[2],  # Orange
        234: accent[4],  # Blue
        237: accent[5],  # Pink
    }
    return {res: manual_map.get(res, 'lightgray') for res in residues_of_interest}


def plot_network(shared_contacts, residues_of_interest, output_prefix, seed=42):
    G = nx.Graph()

    # --- NODE COLORS & OUTLINES ---
    # Catalytic = outline blue, fill white
    # Non-catalytic = outline black, fill white
    node_edgecolors = {}
    node_facecolors = {}

    # --- EDGE COLORS ---
    frst_colors = {
        'NEU': '#555555',   # darker grey
        'MAX': 'red',
        'MIN': 'limegreen'
    }

    # Build graph
    for _, row in shared_contacts.iterrows():
        r1, r2 = row['NumRes1_Ref'], row['NumRes2_Ref']
        l1, l2 = row['MapRes1ID'], row['MapRes2ID']

        # Assign node outline colors
        for r, l in [(r1, l1), (r2, l2)]:
            if r not in G.nodes():
                if r in residues_of_interest:      # catalytic
                    node_edgecolors[r] = "cornflowerblue"
                    node_facecolors[r] = "white"
                else:                              # non-catalytic
                    node_edgecolors[r] = "black"
                    node_facecolors[r] = "white"

                G.add_node(r, label=l)

        G.add_edge(r1, r2, color=frst_colors.get(row['FstConserved'], '#555555'))

    # Circular layout with shuffled nodes
    random.seed(seed)
    nodes = list(G.nodes())
    random.shuffle(nodes)
    angles = np.linspace(0, 2 * np.pi, len(nodes), endpoint=False)
    pos = {node: (np.cos(a) * 2, np.sin(a) * 2) for node, a in zip(nodes, angles)}

    #plt.figure(figsize=(14, 12))
    plt.figure(figsize=(16, 14))

    # --- Draw nodes with outline color ---
    nx.draw_networkx_nodes(
        G, pos,
        node_color=[node_facecolors[n] for n in G.nodes()],
        edgecolors=[node_edgecolors[n] for n in G.nodes()],
        linewidths=6,
        node_size=5000
    )

    # Draw labels
    nx.draw_networkx_labels(
        G, pos,
        labels=nx.get_node_attributes(G, 'label'),
        font_size=20, font_weight='bold'
    )

    # Draw edges
    nx.draw_networkx_edges(
        G, pos,
        edge_color=[G.edges[e]['color'] for e in G.edges()],
        width=6
    )

    # --- Node Legend ---
    node_legend_handles = [
        plt.Line2D([0], [0], marker='o', color='black', label='Non-catalytic',
                   markerfacecolor='white', markeredgewidth=3, markersize=16)
    ]

    node_legend_handles.append(
        plt.Line2D([0], [0], marker='o', color='cornflowerblue',
                   label='Catalytic',
                   markerfacecolor='white', markeredgewidth=3, markersize=16)
    )

    # --- Edge Legend ---
    counts = shared_contacts['FstConserved'].value_counts().to_dict()
    for k in ['MIN', 'NEU', 'MAX']:
        counts.setdefault(k, 0)

    edge_legend_handles = [
        mpatches.Patch(color='red', label=f"MAX ({counts['MAX']})"),
        mpatches.Patch(color='#555555', label=f"NEU ({counts['NEU']})"),
        mpatches.Patch(color='limegreen', label=f"MIN ({counts['MIN']})"),
    ]

    # Add legends
    first_legend = plt.legend(
        handles=node_legend_handles,
        #title='Nodes',
        loc='lower right',
        bbox_to_anchor=(0.99, 0.00),
        fontsize=20,
        frameon=False
    )
    plt.gca().add_artist(first_legend)

    plt.legend(
        handles=edge_legend_handles,
        #title='Frustration states',
        loc='lower right',
        bbox_to_anchor=(0.96, 0.06),
        fontsize=20,
        frameon=False
    )

    plt.axis('off')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_network.svg", dpi=600, bbox_inches='tight')
    plt.close()

    print(f"Network saved as {output_prefix}_network.svg")




def main():
    parser = argparse.ArgumentParser(description="Analyze frustration contacts across 3 files.")
    parser.add_argument("--file1", required=True, help="Path to first file")
    parser.add_argument("--file2", required=True, help="Path to second file")
    parser.add_argument("--file3", required=True, help="Path to third file")
    parser.add_argument("--residues", type=lambda s: {int(item) for item in s.split(",")},
        default={70, 73, 130, 166, 234, 237},
        help="Comma-separated list of residues of interest (e.g., '70,73,130')", required=False)
    parser.add_argument("--prefix", default="output", help="Prefix for output files")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for layout")
    args = parser.parse_args()

    # Step 1: Load and filter
    files = [args.file1, args.file2, args.file3]
    labels = ['Native Family', 'ProteinMPNN 0.1', 'ProteinMPNN 1.0']
    conserved = [read_and_filter(f, args.residues, lbl) for f, lbl in zip(files, labels)]

    # Step 2: Barplot
    plot_barplot(conserved, labels, args.prefix)

    # Step 3: Shared contacts
    sets = [set(df['PairKey']) for df in conserved]
    shared_keys = sets[0] & sets[1] & sets[2]
   
    if not shared_keys:
        print("No conserved contacts shared across all files.")
        return
        
    # New bar plot
    plot_barplot_with_shared(conserved, labels, shared_keys, args.prefix)

    # Merge shared contact rows (just take them from first file where keys match)
    shared_df = conserved[0][conserved[0]['PairKey'].isin(shared_keys)].copy()
    shared_df.to_csv(f"files_contacts_overlay.tsv", sep='\t', index=False)
    
    # Step 4: Network
    plot_network(shared_df, args.residues, args.prefix, seed=args.seed)

if __name__ == "__main__":
    main()
