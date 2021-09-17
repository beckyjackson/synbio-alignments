# 1. BLAST IF against reference genome
# 2. BLAST IF against unicycler assembly

import csv
import logging
import os
import numpy as np
import pandas as pd
import re

from argparse import ArgumentParser
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

BLAST_DB_FILES = [".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"]


def create_blast_db(fasta):
    for ext in BLAST_DB_FILES:
        if not os.path.exists(fasta + ext):
            cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=fasta)
            logging.info(f"creating BLAST database for {fasta}")
            cline()
            break


def get_blast_output(sample_id, genome, genome_fasta, query_fasta):
    seq_id = os.path.basename(os.path.splitext(query_fasta)[0])
    cline = NcbiblastnCommandline(
        query=query_fasta,
        db=genome_fasta,
        evalue=0.001,
        outfmt=5,
        out=f"build/{seq_id}_{genome}.xml",
    )
    logging.info(f"BLASTing {seq_id} against {genome}")
    cline()

    handle = open(f"build/{seq_id}_{genome}.xml")
    record = NCBIXML.read(handle)
    for algn in record.alignments:
        for hsp in algn.hsps:
            return {
                "sample": sample_id,
                "genome": genome,
                "seq_id": seq_id,
                "start": hsp.sbjct_start,
                "end": hsp.sbjct_end,
            }
    return {"sample": sample_id, "genome": genome, "seq_id": seq_id}


def remove_float(i):
    if i is None or pd.isnull(i) or pd.isna(i):
        return None
    return str(i).replace(".0", "")


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "reference_genomes", help="CSV containing master sample ID -> reference genome"
    )
    parser.add_argument("output", help="Path to CSV output")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    else:
        logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    # Read in sample->genome map
    with open(args.reference_genomes, "r") as f:
        reader = csv.reader(f)
        next(reader)
        sample_genome = {row[0]: row[1] for row in reader}

    reference_df = pd.DataFrame()
    assembly_df = pd.DataFrame()
    for sample_id, genome in sample_genome.items():
        sample_dir = f"data/iarpa/TE/{sample_id}"
        seq_dir = os.path.join(sample_dir, "target_sequence")
        if not os.path.exists(seq_dir):
            # logging.info(f"sample {sample_id} contains no target sequences")
            continue
        seq_files = [os.path.join(seq_dir, x) for x in os.listdir(seq_dir)]

        # Find the FASTA for the reference genome
        genome_fasta = None
        for f in os.listdir("data/reference_genomes"):
            if re.match(rf"^{genome}.(fa|fasta|fna)$", f):
                genome_fasta = os.path.join("data/reference_genomes", f)
                break
        if not genome_fasta:
            logging.error(f"unable to find FASTA for {genome} in reference genomes")
            continue

        # Check if we need to create the BLAST database for the ref genome
        create_blast_db(genome_fasta)

        assembly_fasta = os.path.join(sample_dir, "assembly/unicycler_assembly.fasta")
        if not os.path.exists(assembly_fasta) and not os.path.islink(assembly_fasta):
            logging.info(f"no unicycler assembly for sample {sample_id} at " + assembly_fasta)
            continue

        # Check if we need to create the BLAST database for the assembly
        create_blast_db(assembly_fasta)

        for sf in seq_files:
            # BLAST against reference genome
            reference_df = reference_df.append(
                get_blast_output(sample_id, genome, genome_fasta, sf), ignore_index=True
            )
            assembly_df = assembly_df.append(
                get_blast_output(sample_id, "unicycler", assembly_fasta, sf), ignore_index=True
            )
            # TODO: if not aligned to assembly, blast against nt database
            # but we will get lots of 100% matches, how do we filter them?

    output_df = pd.merge(reference_df, assembly_df, on=["sample", "seq_id"], how="outer")
    output_df = output_df.rename(
        columns={
            "genome_x": "ref_genome",
            "start_x": "ref_start",
            "end_x": "ref_end",
            "start_y": "assembly_start",
            "end_y": "assembly_end",
        }
    )
    output_df = output_df[
        [
            "sample",
            "ref_genome",
            "seq_id",
            "ref_start",
            "ref_end",
            "assembly_start",
            "assembly_end",
        ]
    ]
    output_df["ref_start"] = output_df["ref_start"].apply(remove_float)
    output_df["ref_end"] = output_df["ref_end"].apply(remove_float)
    output_df["assembly_start"] = output_df["assembly_start"].apply(remove_float)
    output_df["assembly_end"] = output_df["assembly_end"].apply(remove_float)
    output_df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
