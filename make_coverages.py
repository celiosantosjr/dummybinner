import argparse
import os
import pysam
import pandas as pd


def get_mapped_reads(bam_file):
    """
    Get the total number of mapped reads in a BAM file.
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        return bam.mapped


def create_coverage_table(bam_files, contigs):
    """
    Create a Pandas DataFrame with coverage values for each sample and contig.
    The coverage values are normalized by the total number of mapped reads in each sample.
    """
    # Create an empty DataFrame with columns for each sample
    samples = [os.path.basename(bam_file).replace(".bam", "") for bam_file in bam_files]
    coverage_table = pd.DataFrame(columns=samples, index=contigs)

    # Fill the DataFrame with coverage values
    for bam_file in bam_files:
        sample = os.path.basename(bam_file).replace(".bam", "")
        total_mapped_reads = get_mapped_reads(bam_file)
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for contig in contigs:
                coverage = bam.count_coverage(contig, quality_threshold=0)[0]
                normalized_coverage = [x / total_mapped_reads for x in coverage]
                mean_normalized_coverage = sum(normalized_coverage) / len(normalized_coverage)
                coverage_table.at[contig, sample] = mean_normalized_coverage

    return coverage_table


def write_coverage_table(coverage_table, output_file):
    """
    Write the coverage table to a TSV file.
    """
    coverage_table.to_csv(output_file, sep="\t", header=False, float_format="%.4f", index_label="Contig")


def main():
    parser = argparse.ArgumentParser(description="Create a table of normalized coverage values for each sample and contig.")
    parser.add_argument("bam_files", nargs="+", help="Input BAM files.")
    parser.add_argument("-c", "--contigs", required=True, help="Input FASTA file with reference contigs.")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file with coverage values.")
    args = parser.parse_args()

    # Get a list of contigs from the reference FASTA file
    with pysam.FastaFile(args.contigs) as fasta:
        contigs = fasta.references

    # Create a coverage table with normalized coverage values
    coverage_table = create_coverage_table(args.bam_files, contigs)

    # Write the coverage table to a file
    write_coverage_table(coverage_table, args.output)


if __name__ == "__main__":
    main()

