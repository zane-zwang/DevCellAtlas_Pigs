# Purpose: atlas embedding.
import os
from pathlib import Path

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import to_hex, to_rgb
from matplotlib.lines import Line2D


# 0. Configuration

INPUT = "all_scvi_hvg.h5ad"
OUTDIR = "FigC_output"

os.makedirs(OUTDIR, exist_ok=True)

DPI = 600
FIGSIZE = (8, 8)
POINT_SIZE = 0.6
POINT_ALPHA = 0.85
NUMBER_SIZE = 8

# True: calculate t-SNE
# False: use existing X_tsne
RUN_TSNE = True

plt.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.facecolor": "white"
})


# 1. Load AnnData

adata = sc.read_h5ad(INPUT)

print(adata)
print("Available embeddings:", list(adata.obsm.keys()))

# Correct known spelling error
adata.obs["celltype"] = (
    adata.obs["celltype"]
    .astype("object")
    .replace({"Plasam cell": "Plasma cell"})
    .astype("category")
)


# 2. Calculate t-SNE

if RUN_TSNE:

    if "X_scVI" not in adata.obsm:
        raise KeyError("X_scVI not found in adata.obsm")

    print("\nCalculating t-SNE using X_scVI...")

    sc.tl.tsne(
        adata,
        use_rep="X_scVI",
        perplexity=30,
        random_state=42
    )

    print("t-SNE completed.")

# Check embeddings
for rep in ["X_umap", "X_tsne"]:
    if rep not in adata.obsm:
        raise KeyError(f"{rep} not found.")


# 3. Color dictionaries

ct_cols = {

    # Endothelial
    "Arterial endothelial cell": "#3aa672",
    "Capillary endothelial cell": "#47ad85",
    "Endocardial endothelial cell": "#54b498",
    "Liver endothelial cell": "#61bba1",
    "Lymphatic endothelial cell": "#6ec2b4",
    "Venous endothelial cell": "#7bc9c7",

    # Epithelial
    "BEST4+ epithelial": "#6c3aaa",
    "Brunner's gland cell": "#8B739E",
    "Cholangiocyte": "#7947b5",
    "Ductal-like epithelial cell": "#B88BD8",
    "Enterochromaffin": "#8654c0",
    "Enterocyte": "#9361cb",
    "Enteroendocrine": "#a06ed6",
    "Goblet cell": "#ad7be1",
    "Hepatocyte": "#ba88ec",
    "Intestinal stem cell": "#c795f7",
    "Mesothelial cell": "#d4a2e2",
    "Microfold cell": "#e1afe3",
    "Tuft cell": "#eebce4",

    # Erythroid
    "Erythroblast": "darkslateblue",

    # Immune
    "B cell": "#f53d00",
    "Dendritic cell": "#ff5400",
    "Kupffer cell": "#ff7300",
    "Macrophage": "#fc9a00",
    "Mast cell": "#FFCC00",
    "Microglia": "#FFFF00",
    "Monocyte": "#FFF44F",
    "NK cell": "#FFFF99",
    "NKT cell": "#FFFDD0",
    "Plasma cell": "#fdf7bd",
    "Pre-dendritic cell": "khaki",
    "T cell": "#fde3db",

    # Muscle
    "Cardiomyocyte": "#7f1b1b",
    "Type I myonuclei": "#d96666",
    "Type II a/b myonuclei": "#ec7979",
    "Type II x myonuclei": "#ffb3b3",

    # Neural
    "Astrocyte": "cornflowerblue",
    "Excitatory neuron": "royalblue",
    "Inhibitory neuron": "#1c86ee",
    "Intermediate progenitor cell": "powderblue",
    "Oligodendrocyte": "steelblue",
    "Oligodendrocyte progenitor cell": "dodgerblue",
    "Peripheral glial": "deepskyblue",
    "Peripheral neuron": "lightskyblue",
    "Radial glia": "slateblue",

    # Stromal
    "Adipocyte": "#7f4a1b",
    "Adipose stem and progenitor cell": "#bf6a3f",
    "Adipose-derived mesenchymal stem cell": "peru",
    "Fibro-adipogenic progenitor cell": "#d47f4f",
    "Fibroblast": "#e89a6b",
    "Fibroblast-like cell": "#fcb281",
    "Hepatic stellate cell": "#f0b16f",
    "Interstitial cells of Cajal": "#d59d58",
    "Mesenchymal-like stem cell": "tan",
    "Muscle stem cell": "#7c5a33",
    "Myofibroblast": "#5d4023",
    "Pericyte": "#4c2f0d",
    "Smooth muscle cell": "#bb8e57",
    "Stromal cell": "#9e7532",
    "Tenocyte": "#624510",
    "Vascular leptomeningeal cell": "#785b13",
    "Vascular smooth muscle cell": "goldenrod"
}


# Stage colors
sta_cols = {
    "e55d": "#443a83",
    "e90d": "#31688e",
    "0d": "#21908c",
    "30d": "#35b779",
    "90d": "#8fd744",
    "180d": "#fde725"
}

# Tissue colors
tis_cols = {
    "adipose": "#EEA236FF",
    "cerebrum": "#357EBDFF",
    "duodenum": "#5CB85CFF",
    "heart": "#D43F3A",
    "hypothalamus": "#46B8DAFF",
    "liver": "#20854E99",
    "muscle": "#D43F3AFF"
}

# Lineage colors
lineage_cols = {
    "Endothelial": "#2d9f5f",
    "Epithelial": "#5f2d9f",
    "Erythroid": "#3f3581",
    "Immune": "#fde725",
    "Muscle": "#9f2d2d",
    "Neural": "#2d5f9f",
    "Stromal": "#9f5f2d"
}


# 4. Cell-type grouping

# Explicit lineage assignment
ct_groups = {

    "Endothelial": [
        "Arterial endothelial cell",
        "Capillary endothelial cell",
        "Endocardial endothelial cell",
        "Liver endothelial cell",
        "Lymphatic endothelial cell",
        "Venous endothelial cell"
    ],

    "Epithelial": [
        "BEST4+ epithelial",
        "Brunner's gland cell",
        "Cholangiocyte",
        "Ductal-like epithelial cell",
        "Enterochromaffin",
        "Enterocyte",
        "Enteroendocrine",
        "Goblet cell",
        "Hepatocyte",
        "Intestinal stem cell",
        "Mesothelial cell",
        "Microfold cell",
        "Tuft cell"
    ],

    "Erythroid": [
        "Erythroblast"
    ],

    "Immune": [
        "B cell",
        "Dendritic cell",
        "Kupffer cell",
        "Macrophage",
        "Mast cell",
        "Microglia",
        "Monocyte",
        "NK cell",
        "NKT cell",
        "Plasma cell",
        "Pre-dendritic cell",
        "T cell"
    ],

    "Muscle": [
        "Cardiomyocyte",
        "Type I myonuclei",
        "Type II a/b myonuclei",
        "Type II x myonuclei"
    ],

    "Neural": [
        "Astrocyte",
        "Excitatory neuron",
        "Inhibitory neuron",
        "Intermediate progenitor cell",
        "Oligodendrocyte",
        "Oligodendrocyte progenitor cell",
        "Peripheral glial",
        "Peripheral neuron",
        "Radial glia"
    ],

    "Stromal": [
        "Adipocyte",
        "Adipose stem and progenitor cell",
        "Adipose-derived mesenchymal stem cell",
        "Fibro-adipogenic progenitor cell",
        "Fibroblast",
        "Fibroblast-like cell",
        "Hepatic stellate cell",
        "Interstitial cells of Cajal",
        "Mesenchymal-like stem cell",
        "Muscle stem cell",
        "Myofibroblast",
        "Pericyte",
        "Smooth muscle cell",
        "Stromal cell",
        "Tenocyte",
        "Vascular leptomeningeal cell",
        "Vascular smooth muscle cell"
    ]
}


# 5. Sort and validate

lineage_order = list(ct_groups.keys())

# Keep only cell types present in the actual AnnData
actual_ct = set(
    adata.obs["celltype"].dropna().astype(str).unique()
)

# Sort alphabetically within each lineage
for lineage in lineage_order:

    ct_groups[lineage] = sorted(
        [
            ct for ct in ct_groups[lineage]
            if ct in actual_ct
        ],
        key=str.casefold
    )

# Flatten the grouped cell types
ordered_ct = [
    ct
    for lineage in lineage_order
    for ct in ct_groups[lineage]
]

# Check completeness
missing = actual_ct - set(ordered_ct)

if missing:
    raise ValueError(
        f"Cell types missing from ct_groups: {sorted(missing)}"
    )

if len(ordered_ct) != len(set(ordered_ct)):
    raise ValueError("Duplicate cell types in ct_groups.")

# Check color coverage
missing_colors = set(ordered_ct) - set(ct_cols)

if missing_colors:
    raise ValueError(
        f"Missing cell-type colors: {sorted(missing_colors)}"
    )

# Check tissue and stage colors
for col, palette in [
    ("tissue", tis_cols),
    ("stage", sta_cols)
]:
    values = set(
        adata.obs[col].dropna().astype(str).unique()
    )

    missing_values = values - set(palette)

    if missing_values:
        raise ValueError(
            f"Missing {col} colors: {sorted(missing_values)}"
        )

# Generate continuous numbering
ct_number = {
    ct: i + 1
    for i, ct in enumerate(ordered_ct)
}

# Export legend table
legend_records = []

for lineage in lineage_order:

    for ct in ct_groups[lineage]:

        legend_records.append({
            "Number": ct_number[ct],
            "Celltype": ct,
            "Celllineage": lineage,
            "Color": to_hex(ct_cols[ct])
        })

legend_df = pd.DataFrame(legend_records)

legend_df.to_csv(
    f"{OUTDIR}/FigC_celltype_legend.csv",
    index=False
)

# Set Scanpy colors
adata.obs["celltype"] = pd.Categorical(
    adata.obs["celltype"],
    categories=ordered_ct,
    ordered=True
)

adata.uns["celltype_colors"] = [
    to_hex(ct_cols[ct])
    for ct in ordered_ct
]

print(f"\nTotal cell types: {len(ordered_ct)}")
print("Color validation passed.")


# 6. Common figure saving function

def save_figure(fig, name, png=True):

    base = os.path.join(OUTDIR, name)

    if png:
        fig.savefig(
            base + ".png",
            dpi=DPI,
            facecolor="white"
        )

    fig.savefig(
        base + ".pdf",
        dpi=DPI,
        facecolor="white"
    )

    fig.savefig(
        base + ".svg",
        dpi=DPI,
        facecolor="white"
    )

    plt.close(fig)


# 7. Common embedding plotting function

def plot_embedding(
    embedding="umap",
    color_by="celltype",
    numbered=False
):

    coords = np.asarray(
        adata.obsm[f"X_{embedding}"]
    )[:, :2]

    labels = (
        adata.obs[color_by]
        .astype("object")
        .to_numpy()
    )

    if not np.isfinite(coords).all():
        raise ValueError(
            f"Non-finite coordinates in X_{embedding}"
        )

    # Select palette and plotting order
    if color_by == "celltype":
        palette = ct_cols
        categories = ordered_ct

    elif color_by == "tissue":
        palette = tis_cols
        categories = list(tis_cols)

    elif color_by == "stage":
        palette = sta_cols
        categories = list(sta_cols)

    else:
        raise ValueError(f"Unknown color_by: {color_by}")

    # Create a square figure
    fig, ax = plt.subplots(
        figsize=FIGSIZE,
        dpi=150
    )

    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Draw cells
    for category in categories:

        idx = labels == category

        if not np.any(idx):
            continue

        xy = coords[idx]

        ax.scatter(
            xy[:, 0],
            xy[:, 1],
            s=POINT_SIZE,
            c=[palette[category]],
            alpha=POINT_ALPHA,
            linewidths=0,
            rasterized=True,
            zorder=1
        )

    # Add numbered cell-type labels

    if numbered:

        for ct in ordered_ct:

            xy = coords[labels == ct]

            if len(xy) == 0:
                continue

            # Median coordinate
            center = np.median(xy, axis=0)

            # Find nearest actual cell
            distance = np.sum(
                (xy - center) ** 2,
                axis=1
            )

            x, y = xy[np.argmin(distance)]

            # Number with white outline
            ax.text(
                x,
                y,
                str(ct_number[ct]),
                fontsize=NUMBER_SIZE,
                fontweight="bold",
                color="black",
                ha="center",
                va="center",
                zorder=10,
                bbox=dict(
                    boxstyle="circle,pad=0.15",
                    facecolor=ct_cols[ct],
                    edgecolor="white",
                    linewidth=0.6,
                    alpha=0.98
                )
            )

    # Make plotting area square (1:1)

    xmin, ymin = coords.min(axis=0)
    xmax, ymax = coords.max(axis=0)

    cx = (xmin + xmax) / 2
    cy = (ymin + ymax) / 2

    span = max(
        xmax - xmin,
        ymax - ymin
    ) * 1.06

    ax.set_xlim(
        cx - span / 2,
        cx + span / 2
    )

    ax.set_ylim(
        cy - span / 2,
        cy + span / 2
    )

    ax.set_aspect("equal")
    ax.set_box_aspect(1)

    # No axis, no legend
    ax.axis("off")

    # Fill figure while preserving square aspect
    fig.subplots_adjust(
        left=0.01,
        right=0.99,
        bottom=0.01,
        top=0.99
    )

    # Save

    if numbered:
        name = f"FigC_{embedding}_numbered"
    else:
        name = f"FigC_{embedding}_{color_by}"

    save_figure(fig, name)

    print("Saved:", name)


# 8. Independent cell-type color legend

def plot_celltype_legend():

    # Three-column layout:
    # 1: Endothelial + Epithelial + Erythroid + Immune
    # 2: Muscle + Neural
    # 3: Stromal

    columns = [
        [
            "Endothelial",
            "Epithelial",
            "Erythroid",
            "Immune"
        ],
        [
            "Muscle",
            "Neural"
        ],
        [
            "Stromal"
        ]
    ]

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(17, 9)
    )

    # Determine the required number of rows per column
    def n_rows(groups):
        return sum(
            2 + len(ct_groups[g])
            for g in groups
        )

    max_rows = max(
        n_rows(groups)
        for groups in columns
    )

    for ax, groups in zip(axes, columns):

        ax.set_xlim(0, 1)
        ax.set_ylim(max_rows + 1, -1)
        ax.axis("off")

        y = 0

        for lineage in groups:

            names = ct_groups[lineage]

            if not names:
                continue

            # Lineage header
            ax.text(
                0.02,
                y,
                lineage,
                fontsize=12,
                fontweight="bold",
                color=lineage_cols[lineage],
                ha="left",
                va="center"
            )

            y += 1.3

            # Cell types
            for ct in names:

                number = ct_number[ct]

                # Colored circle
                ax.scatter(
                    0.055,
                    y,
                    s=170,
                    color=ct_cols[ct],
                    edgecolor="none",
                    zorder=2
                )

                # Number inside circle
                ax.text(
                    0.055,
                    y,
                    str(number),
                    fontsize=6.5,
                    fontweight="bold",
                    color="black",
                    ha="center",
                    va="center",
                    zorder=3
                )

                # Cell-type name
                ax.text(
                    0.105,
                    y,
                    ct,
                    fontsize=9,
                    color="#222222",
                    ha="left",
                    va="center"
                )

                y += 1

            y += 0.8

    fig.subplots_adjust(
        left=0.02,
        right=0.99,
        top=0.98,
        bottom=0.02,
        wspace=0.12
    )

    save_figure(
        fig,
        "FigC_cell_color_legend",
        png=False
    )

    print("Saved: FigC_cell_color_legend")


# 9. Generate all figures

print("\nGenerating Figure C...")

# Cell-type numbered plots
plot_embedding(
    embedding="umap",
    color_by="celltype",
    numbered=True
)

plot_embedding(
    embedding="tsne",
    color_by="celltype",
    numbered=True
)

# Tissue-colored plots
plot_embedding(
    embedding="umap",
    color_by="tissue"
)

plot_embedding(
    embedding="tsne",
    color_by="tissue"
)

# Stage-colored plots
plot_embedding(
    embedding="umap",
    color_by="stage"
)

plot_embedding(
    embedding="tsne",
    color_by="stage"
)

# Independent legend
plot_celltype_legend()

print("\nAll figures completed.")
print("Output directory:", OUTDIR)
