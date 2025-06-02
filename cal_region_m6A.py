import polars as pl
import re, argparse, os, subprocess

parser = argparse.ArgumentParser(description="Calculate region/gene level m6A.")

parser.add_argument("--input", "-i", nargs="+", help="unfiltered *totalm6A.FDR.csv files")
parser.add_argument("--Afile", "-A", 
                  help="*m6A.used file from merge_multisam_m6A.py, which recording used A sites to count as denominator")
parser.add_argument("--genome", "-g", default="hg38", help="genome used, (hg38 [default] or mm10")
parser.add_argument("--output", "-o", help="output file path")

args = parser.parse_args()
fdist = args.Afile + "_dist.measures.txt"

# annotate used A file
if not os.path.exists(fdist):
  subprocess.call("annotate_sites.sh -s " + args.Afile + " -m m6A -g " + args.genome, shell=True)


# read dist info for used sites
df_dist = pl.read_csv(fdist, separator = "\t")

# get info about whether sites passed in any sample
df_merged = pl.read_csv(args.Afile, separator = "\t")
df_merged.columns = ["chr", "coord", "Strand", "Gene", "Passed"]
df_dist = df_dist.join(
    df_merged, on = ["chr", "coord"], how = "left", coalesce = True
).with_columns(
    pl.when(pl.col("rel_location") <= 1)
    .then(pl.lit("utr5"))
    .when(pl.col("rel_location") <= 2)
    .then(pl.lit("cds"))
    .otherwise(pl.lit("utr3"))
    .alias("Region")
)
# df_dist


for fin in args.input:
    sam = re.match("(.*/)?([^/]+).totalm6A.FDR.csv$", fin).group(2)
    df_sam = pl.read_csv(fin, separator = "\t")

    # merge with dist info
    df = df_dist.join(
        df_sam.rename({"Chr": "chr", "Sites": "coord"}),
        on = ["chr", "coord", "Strand", "Gene"], 
        how = "left", coalesce = True
    ).select(
        pl.col("chr", "coord", "Strand", "gene_name", "Region", "refseqID", "Passed"),
        pl.col("AGcov").alias("AGcov_"+sam),
        pl.col("Acov").alias("Acov_"+sam),
    ).fill_null(False)
    df

    # calculate region level m6A ratio
    df_region = df.group_by(["gene_name", "Region"]).agg(
        pl.col("AGcov_"+sam).sum(),
        pl.col("Acov_"+sam).filter(pl.col("Passed")).sum()
    ).with_columns(
        (pl.col("Acov_"+sam)/pl.col("AGcov_"+sam)*100).alias("Ratio_"+sam)
    )

    # calculate gene level m6A ratio
    df_gene = df.group_by("gene_name").agg(
        pl.lit("gene").alias("Region"),
        pl.col("AGcov_"+sam).sum(),
        pl.col("Acov_"+sam).filter(pl.col("Passed")).sum()
    ).with_columns(
        (pl.col("Acov_"+sam)/pl.col("AGcov_"+sam)*100).alias("Ratio_"+sam)
    )
    df_region = pl.concat([df_region, df_gene]).sort(["gene_name", "Region"])

    # merge all samples
    if fin == args.input[0]:
        df_final = df_region
    else:
        df_final = df_final.join(
            df_region, 
            on=["gene_name", "Region"], how="inner", coalesce=True
        )

# output merged list
df_final.write_csv(args.output, separator='\t')