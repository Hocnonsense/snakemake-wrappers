"""Snakemake wrapper for [metapop](https://github.com/metaGmetapop/metapop)."""

__author__ = "Aoran Hu"
__copyright__ = "Copyright 2025, Aoran Hu"
__email__ = "hwrn.aou@sjtu.edu.cn"
__license__ = "MIT"


from pathlib import Path
from tempfile import TemporaryDirectory
from snakemake.shell import shell
from snakemake.script import snakemake

import pandas as pd


def mkparent(path):
    p = Path(path)
    if not p.parent.exists():
        p.parent.mkdir(parents=True)
    elif p.exists():
        shell(f"rm -r {p}")


def either_output(name):
    if name in snakemake.output.keys():
        yield snakemake.output[name]
    if name in snakemake.params.keys():
        yield snakemake.params[name]
    if name in snakemake.log.keys():
        yield snakemake.log[name]


folder_to_delete: list[TemporaryDirectory] = []
genome_file_suffix = {".fasta", ".fas", ".fa", ".fna", ".fsa_nt"}
if "suffix" in snakemake.params.keys():
    genome_file_suffix = set(genome_file_suffix)


def make_input_folder(folder_or_lsfile: Path):
    """
    check input file or folder, and return a folder to fit software input.

    accept inputs:
    1. collect input genomes from:
        1.1. a folder -> all genomes in the folder with given suffix;
        1.2. a file that list genomes:
            2.1 plain list of files -> try to infer genome name from the file names;
            2.2 a tsv with two columns witout header -> extract name[tab]path pairs;

            use "#" at the start of a line to comment out the line.
            will check if all files exist

    2. create a temp folder with soft links to the files.

    For softlinking, the file with duplicated basenames will cause problems.

    return: folder_path, suffix
        As a temp folder, the folder should be deleted after use (stored in `folder_to_delete`).
    """
    extra_suffix = list(genome_file_suffix)[0] if len(genome_file_suffix) == 1 else ""
    if folder_or_lsfile.is_dir():
        files_dict = {
            i.name: i
            for i in folder_or_lsfile.glob("*")
            if i.is_file() and any(i.name.endswith(j) for j in genome_file_suffix)
        }
    else:
        with open(folder_or_lsfile) as f:
            files = {i.strip() for i in f if not i.startswith("#")}
            if all("\t" in str(i) for i in files):
                extra_suffix = extra_suffix or ".fa"
                files_dict = {
                    k: Path(f"{v}{extra_suffix}")
                    for k, v in (j.split("\t")[:2] for j in files)
                }
            else:
                files_dict = {i.name: i for i in (Path(j) for j in files)}
        not_exist_files = [str(i) for i in files_dict.values() if not i.is_file()]
        if not_exist_files:
            raise ValueError(f"Files does not exist: \n" + "\n".join(not_exist_files))
    tempdir = TemporaryDirectory(delete=False)
    # fine the basename and check if there are duplicated basenames
    for name, path in files_dict.items():
        shell(f"ln -rs {path} {tempdir.name}/{name}{extra_suffix}")
    folder_to_delete.append(tempdir)
    return tempdir.name, extra_suffix


run_temp = TemporaryDirectory(delete=False)
folder_to_delete.append(run_temp)

genomes_dir, extra_suffix = make_input_folder(Path(snakemake.input.genomes))
pyanidb = f"{run_temp.name}/pyanidb"
pyani_exec = snakemake.params.methods

shell(
    f"""
    mkdir -p {run_temp.name}/report

    pyani createdb --dbpath {pyanidb}
    pyani index -i {genomes_dir}
    pyani {pyani_exec} --dbpath {pyanidb} -i {genomes_dir} -o {run_temp.name} \\
        --workers {snakemake.threads}

    pyani report --dbpath {pyanidb} -o {run_temp.name}/report \\
        --genomes \\
        --run_matrices 1
    """,
)


def read_pyani_out(output_base: Path, runid=1):
    gid = pd.read_csv(
        output_base / "report" / "genomes.tab",
        sep="\t",
        index_col=0,
    )
    id2name = gid.assign(
        index=lambda df: "Genome_id:" + df["genome ID"].astype(str),
        name=lambda df: df.path.str.rsplit("/", n=1).apply(lambda x: x[-1]),
    ).set_index("index")["name"]
    if extra_suffix:
        id2name = id2name.apply(lambda x: x[: -len(extra_suffix)])

    df: pd.DataFrame = None
    for file in (output_base / "report").glob(f"matrix_*_{runid}.tab"):
        df1 = (
            pd.read_csv(file, sep="\t", index_col=0)
            .reset_index()
            .melt(id_vars=["index"], var_name="genome_2", value_name=file.name[7:-6])
        )
        if df is None:
            df = df1
        else:
            df = df.merge(df1, on=["index", "genome_2"])
    df_out = df.assign(
        index=df["index"].map(id2name),
        genome_2=df["genome_2"].map(id2name),
    ).rename(columns={"index": "genome_1"})
    return df_out


for pyani_tsv in either_output("tsv"):
    mkparent(pyani_tsv)
    df_out = read_pyani_out(Path(run_temp.name), 1)
    df_out.to_csv(pyani_tsv, sep="\t", index=False)

for pyanidb_out in either_output("pyanidb"):
    mkparent(pyanidb_out)
    shell(f"mv {pyanidb} {pyanidb_out}")

for pyani_report in either_output("report"):
    mkparent(pyani_report)
    shell(f"mv {run_temp.name}/report {pyani_report}")

for pyani_outdir in either_output("outdir"):
    mkparent(pyani_outdir)
    shell(f"mv {run_temp.name} {pyani_outdir}")


for folder in folder_to_delete:
    shell(f"rm -rf {folder.name}")
