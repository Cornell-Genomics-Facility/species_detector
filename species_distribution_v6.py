"""
This script processes BLAST output files in the specified directory to analyze species distributions
and detect significant differences in distributions across files.

Usage example:
    python species_distribution_v5.py --base_directory /path/to/directory --ecut 1e-20

Author: Paul Munn, Genomics Innovation Hub, Cornell University

Version history:
    - 08/01/2024: Original version
    - 08/05/2024: Create html file from stacked bar plot
    - 08/13/2024: Add a stacked bar plot to group each species by order (e.g. "Primates", "Carnivora", "Rodentia", "Chiroptera", etc.)
    - 09/05/2024: Add hover text to stacked bar plots
"""

import os
import argparse
import pandas as pd
from collections import defaultdict
import json
from Bio import Entrez
import plotly.graph_objects as go
import re

# Entrez email configuration
Entrez.email = "prm88@cornell.edu"  # Update with your email address

# Known Phyla and Orders
known_phylums = ["Anthocerotophyta", "Bryophyta", "Chlorophyta", "Charophyta", "Cycadophyta", "Ginkgophyta", "Gnetophyta", "Hepatophyta", "Lycopodiophyta", "Magnoliophyta", "Marchantiophyta", "Pinophyta", "Pteridophyta", "Porifera", "Cnidaria", "Ctenophora", "Platyhelminthes", "Nematoda", "Mollusca", "Annelida", "Arthropoda", "Echinodermata", "Chordata", "Nemertea", "Bryozoa", "Rotifera", "Proteobacteria", "Firmicutes", "Actinobacteria", "Bacteroidetes", "Cyanobacteria", "Spirochaetes", "Chlamydiae", "Planctomycetes", "Verrucomicrobia", "Acidobacteria", "Aquificae", "Thermotogae", "Chlorobi", "Chloroflexi"]
known_orders = ["Alismatales", "Apiales", "Asterales", "Brassicales", "Caryophyllales", "Cornales", "Ericales", "Fabales", "Fagales", "Lamiales", "Malpighiales", "Malvales", "Myrtales", "Poales", "Rosales", "Sapindales", "Solanales", "Primates", "Carnivora", "Cetacea", "Chiroptera", "Rodentia", "Proboscidea", "Lagomorpha", "Perissodactyla", "Artiodactyla", "Sirenia", "Diptera", "Coleoptera", "Lepidoptera", "Hymenoptera", "Pseudomonadales", "Enterobacteriales", "Lactobacillales", "Clostridiales", "Actinomycetales", "Bacteroidales", "Nostocales", "Chroococcales", "Spirochaetales", "Chlamydiales", "Planctomycetales", "Verrucomicrobiales", "Acidobacteriales"]


# Function to fetch lineage information from NCBI Entrez
def fetch_lineage(taxid):
    return "Unknown", "Unknown"  # Placeholder for future implementation
    # """Fetch the lineage of a taxon ID using NCBI Entrez."""
    # handle = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
    # records = Entrez.read(handle)
    # handle.close()
    # lineage = records[0]["LineageEx"]
    # phylum = None
    # order = None
    # for taxon in lineage:
    #     if taxon["Rank"] == "phylum" and taxon["ScientificName"] in KNOWN_PHYLUMS:
    #         phylum = taxon["ScientificName"]
    #     if taxon["Rank"] == "order" and taxon["ScientificName"] in KNOWN_ORDERS:
    #         order = taxon["ScientificName"]
    # return phylum, order


def parse_blast_file(file_path, ecut, seen_identifiers):
    """Parse a BLAST file and collect species counts based on a given significance cutoff."""
    species_count = defaultdict(int)
    blastname_count = defaultdict(int)
    total_count = 0
    taxonomy_ids = {}
    kingdom = "Unknown"

    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            read_id = columns[0]
            taxonomy_id = columns[1]
            significance = float(columns[10])
            species = columns[13]
            blastname = columns[15]
            kingdom = columns[16]

            # Ignore rows with seen identifiers or significance above the cutoff
            # print('Read ID:', read_id)
            # print('Species:', species)
            # print('BlastName:', blastname)
            # print('Significance:', significance)
            if read_id in seen_identifiers or significance >= ecut:
                continue

            seen_identifiers.add(read_id)
            species_count[species] += 1
            blastname_count[blastname] += 1
            total_count += 1
            taxonomy_ids[species] = taxonomy_id  # Store taxonomy id for lineage info

    return species_count, blastname_count, total_count, taxonomy_ids, kingdom


def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split('(\d+)', s)]


def generate_output(base_directory, blast_directory, ecut, sort_by):
    """Main function to orchestrate parsing and generating outputs."""
    blast_files = [f for f in os.listdir(base_directory + '/' + blast_directory) if f.endswith('.blast')]
    if sort_by.lower() == "name":
        blast_files.sort(key=natural_sort_key, reverse=True)  # Sort files by name
    all_species_counts = defaultdict(int)
    all_blastname_counts = defaultdict(int)
    file_species_counts = []
    file_blastname_counts = []

    # Parse each BLAST file and aggregate species and blastname counts
    for blast_file in blast_files:
        # print('Blast File:', blast_file)
        file_path = os.path.join(base_directory, blast_directory, blast_file)
        seen_identifiers = set()
        species_count, blastname_count, total_count, taxonomy_ids, kingdom = parse_blast_file(file_path, ecut,
                                                                                              seen_identifiers)

        file_species_counts.append({
            "filename": blast_file,
            "species_count": species_count,
            "blastname_count": blastname_count,
            "total_count": total_count,
            "taxonomy_ids": taxonomy_ids,
            "kingdom": kingdom
        })
        # print('Species Count:', species_count)
        # print('BlastName Count:', blastname_count)
        # print('Total Count:', total_count)
        # print('Taxonomy IDs:', taxonomy_ids)
        # print('Kingdom:', kingdom)

        for species, count in species_count.items():
            all_species_counts[species] += count
        for blastname, count in blastname_count.items():
            all_blastname_counts[blastname] += count

    total_species_sum = sum(all_species_counts.values())
    species_avg_distribution = {species: count / total_species_sum for species, count in all_species_counts.items()}

    total_blastname_sum = sum(all_blastname_counts.values())
    blastname_avg_distribution = {blastname: count / total_blastname_sum for blastname, count in
                                  all_blastname_counts.items()}

    # Write output file
    with open(os.path.join(base_directory, "species_stats.txt"), 'w') as out_file:
        out_file.write("\t".join(["Filename", "Blastname", "TaxonomyID", "Kingdom", "Species", "TotalCount",
                                  "SpeciesPercentage", "TotalPercentageAcrossFiles", "AverageDistribution"]) + "\n")

        for file_stats in file_species_counts:
            blast_file = file_stats["filename"]
            species_count = file_stats["species_count"]
            total_count = file_stats["total_count"]
            taxonomy_ids = file_stats["taxonomy_ids"]
            kingdom = file_stats["kingdom"]

            for species, count in species_count.items():
                species_percentage = count / total_count
                total_percentage_across_files = all_species_counts[species] / total_species_sum
                avg_distribution = species_avg_distribution[species]

                taxonomy_id = taxonomy_ids[species]

                out_file.write("\t".join(
                    [blast_file, json.dumps(file_stats["blastname_count"]), taxonomy_id, kingdom, species, str(count),
                     f"{species_percentage:.4f}", f"{total_percentage_across_files:.4f}",
                     f"{avg_distribution:.4f}"]) + "\n")

    # Generate HTML Report
    html_report = "<html><body>"
    overall_top_species = defaultdict(int)
    overall_top_blastname = defaultdict(int)

    for file_stats in file_species_counts:
        blast_file = file_stats["filename"]
        species_count = file_stats["species_count"]
        sorted_species = sorted(species_count.items(), key=lambda item: item[1], reverse=True)[:5]
        for species, count in sorted_species:
            overall_top_species[species] += count

        blastname_count = file_stats["blastname_count"]
        sorted_blastname = sorted(blastname_count.items(), key=lambda item: item[1], reverse=True)[:5]
        for blastname, count in sorted_blastname:
            overall_top_blastname[blastname] += count

    sorted_overall_species = sorted(overall_top_species.items(), key=lambda item: item[1], reverse=True)[:5]
    sorted_overall_blastname = sorted(overall_top_blastname.items(), key=lambda item: item[1], reverse=True)[:5]

    print("Overall Top Species:", sorted_overall_species[0])
    print("Overall Top Blastname:", sorted_overall_blastname[0])
    # print('File species counts:', file_species_counts[:5])

    if sort_by.lower() == "topspecies":
        # If not sorted by name, sort by the top species count
        species_to_sort_by = sorted_overall_species[0][0]
        file_species_counts = sorted(
            file_species_counts,
            key=lambda x: x['species_count'].get(species_to_sort_by, 0) / x['total_count'],
            reverse=False)  # Use reverse=True if you want descending order

    html_report += "<h1>Overall Top 5 Species</h1><ul>"
    for species, count in sorted_overall_species:
        html_report += f"<li>{species}: {count} ({count / total_species_sum:.2%})</li>"
    html_report += "</ul>"

    html_report += "<h1>Top 5 Species per BLAST file</h1>"
    for file_stats in file_species_counts:
        blast_file = file_stats["filename"]
        species_count = file_stats["species_count"]
        sorted_species = sorted(species_count.items(), key=lambda item: item[1], reverse=True)[:5]

        html_report += f"<h2>{blast_file}</h2><ul>"
        for species, count in sorted_species:
            html_report += f"<li>{species}: {count} ({count / file_stats['total_count']:.2%})</li>"
        html_report += "</ul>"

    html_report += "</body></html>"
    with open(os.path.join(base_directory, "species_summary_mqc.html"), 'w') as out_file:
        out_file.write(html_report)

    # Generate species distribution plot
    bar_height = 25  # Height of each bar in pixels
    plot_height = max(1000, bar_height * len(file_species_counts))
    figure_data = []
    for file_stats in file_species_counts:
        blast_file = file_stats["filename"]
        species_count = file_stats["species_count"]
        for species, count in species_count.items():
            species_percentage = count / file_stats["total_count"]
            if species_percentage > 0.001:  # Only include if percentage > 0.1%
                figure_data.append((blast_file, species, species_percentage))

    df = pd.DataFrame(figure_data, columns=["File", "Species", "Percentage"])

    fig1 = go.Figure()
    species_order = df.groupby("Species")["Percentage"].sum().sort_values(ascending=False).index.tolist()
    # First, get all unique combinations of files and species
    all_files = df["File"].unique()
    all_species = species_order  # We already have this from your existing code

    # Create a complete DataFrame with all combinations
    index = pd.MultiIndex.from_product([all_files, all_species], names=["File", "Species"])
    df_complete = pd.DataFrame(index=index).reset_index()

    # Merge with original data, filling NaN with 0
    df_complete = df_complete.merge(df, on=["File", "Species"], how="left")
    df_complete["Percentage"] = df_complete["Percentage"].fillna(0)

    # Create the stacked bar plot
    for species in species_order:
        species_data = df_complete[df_complete["Species"] == species]
        fig1.add_trace(go.Bar(
            y=species_data["File"].str.rsplit('.', n=1).str[0],
            x=species_data["Percentage"],
            name=f"{species[:45]} ({round(all_species_counts[species] * 100 / total_species_sum, 2)}%)",
            orientation='h',
            text=species,
            hoverinfo='x+y+text',
        ))

    fig1.update_layout(
        barmode='stack',
        title='Species Distribution by Sample',
        xaxis_title='Percentage',
        yaxis_title='Sample',
        xaxis=dict(tickformat=".2%"),
        bargap=0.2,
        bargroupgap=0.0,
        height=plot_height,
        width=1400,
        legend=dict(traceorder='normal', title="Species (Overall Percentage)", itemsizing='constant', orientation='v')
    )

    fig1.write_html(os.path.join(base_directory, "species_stacked_barplot_mqc.html"))
    fig1.show()

    if sort_by.lower() == "topspecies":
        # If not sorted by name, sort by the top blastname count
        blastname_to_sort_by = sorted_overall_blastname[0][0]
        file_species_counts = sorted(
            file_species_counts,
            key=lambda x: x['blastname_count'].get(blastname_to_sort_by, 0) / x['total_count'],
            reverse=False)  # Use reverse=True if you want descending order

    # Generate blastname distribution plot
    figure_data = []
    for file_stats in file_species_counts:
        blast_file = file_stats["filename"]
        blastname_count = file_stats["blastname_count"]
        for blastname, count in blastname_count.items():
            blastname_percentage = count / file_stats["total_count"]
            if blastname_percentage > 0.001:  # Only include if percentage > 0.1%
                figure_data.append((blast_file, blastname, blastname_percentage))

    df2 = pd.DataFrame(figure_data, columns=["File", "Blastname", "Percentage"])

    fig2 = go.Figure()
    blastname_order = df2.groupby("Blastname")["Percentage"].sum().sort_values(ascending=False).index.tolist()
    # Do the same for blastnames
    all_blastnames = blastname_order  # We already have this from your existing code

    # Create a complete DataFrame with all combinations
    index = pd.MultiIndex.from_product([all_files, all_blastnames], names=["File", "Blastname"])
    df2_complete = pd.DataFrame(index=index).reset_index()

    # Merge with original data, filling NaN with 0
    df2_complete = df2_complete.merge(df2, on=["File", "Blastname"], how="left")
    df2_complete["Percentage"] = df2_complete["Percentage"].fillna(0)

    for blastname in blastname_order:
        blastname_data = df2_complete[df2_complete["Blastname"] == blastname]
        fig2.add_trace(go.Bar(
            y=blastname_data["File"].str.rsplit('.', n=1).str[0],
            x=blastname_data["Percentage"],
            name=f"{blastname[:45]} ({round(all_blastname_counts[blastname] * 100 / total_blastname_sum, 2)}%)",
            orientation='h',
            text=blastname,
            hoverinfo='x+y+text',
        ))

    fig2.update_layout(
        barmode='stack',
        title='Blastname Distribution by Sample',
        xaxis_title='Percentage',
        yaxis_title='Sample',
        xaxis=dict(tickformat=".2%"),
        bargap=0.2,
        bargroupgap=0.0,
        height=plot_height,
        width=1300,
        legend=dict(traceorder='normal', title="Blastname (Overall Percentage)", itemsizing='constant', orientation='v')
    )

    fig2.write_html(os.path.join(base_directory, "blastname_stacked_barplot_mqc.html"))
    fig2.show()


def main():
    parser = argparse.ArgumentParser(description="Process BLAST search results for species distribution analysis.")
    parser.add_argument('--base_directory', type=str, help='Path to the working directory containing BLAST files.',
                        default='C:/Users/prm88/Documents/Box/genome_innovation_hub/genomics_facility/fastq_species_detector/Project_14242_32832/Species_detector_report')
                        # default='Species_detector_report')
    parser.add_argument('--blast_directory', type=str, default='blast_out_dir', help='Path to the working directory containing BLAST files.')
    parser.add_argument('--ecut', type=float, default=1e-5, help='E-value cutoff threshold for filtering BLAST results.')
    parser.add_argument('--sort_by', type=str, help="Sorting for samples", default="name")
    args = parser.parse_args()

    print(f"Base Directory: {args.base_directory}")
    print(f"Blast Directory: {args.blast_directory}")
    print(f"E-value Cutoff: {args.ecut}")

    generate_output(args.base_directory, args.blast_directory, args.ecut, args.sort_by)


if __name__ == "__main__":
    main()
