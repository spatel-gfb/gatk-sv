#!/bin/python

# Synopsis:
#  Creates input values for a new test batch from a GATKSVPipelineBatch run
#

import argparse
import json
import sys


INPUT_KEYS = set([
    "name",
    "samples",
    "bam_or_cram_files",
    "requester_pays_crams",
    "gvcfs",
    "ped_file",
    "contig_ploidy_model_tar",
    "gcnv_model_tars",
    "qc_definitions",
    "outlier_cutoff_table"
])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata", help="GATKSVPipelineBatch metadata JSON file")
    args = parser.parse_args()

    with open(args.metadata, 'r') as f:
        metadata = json.load(f)
    values = {key.replace("GATKSVPipelineBatch.", ""): value for key, value in metadata["outputs"].items()
              if value is not None}
    inputs = metadata["inputs"]
    for raw_key in set(inputs.keys()).intersection(INPUT_KEYS):
        key = raw_key.split('.')[-1]
        values[key] = inputs[key]
    for key in INPUT_KEYS - set(values.keys()):
        sys.stderr.write(f"Warning: expected workflow input '{key}' not found in metadata. You will need to add "
                         f"this entry manually.\n")
        values[key] = None

    print(json.dumps(values, sort_keys=True, indent=4))


if __name__ == "__main__":
    main()
