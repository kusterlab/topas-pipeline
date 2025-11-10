# compute substrate phosphorylation scores (a.k.a. kinase activity) for cytoplasmic
# kinases using the confident substrates from Florian Bayer's decryptM experiments

import sys
import logging
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np
from tqdm import tqdm

import psite_annotation as pa

from .. import sample_metadata

tqdm.pandas()

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)

ATR_seeds = pd.DataFrame(
    {
        "Modified sequence": sorted(
            set(
                [
                    "HEAISS(ph)QESK",
                    "LIDIVSS(ph)QK",
                    "YSSS(ph)QPEPR",
                    "ETPS(ph)QENGPTAK",
                    "KQPPVSPGTALVGS(ph)QK",
                    "SS(ph)QPLASK",
                    "IAGMS(ph)QK",
                    "SSGISS(ph)QNSSTSDGDR",
                    "SVSS(ph)QSSSSVSSQVTTAGSGK",
                    "GSHIS(ph)QGNEAEER",
                    "ELGTGLS(ph)QKR",
                    "DALAALETPGRPS(ph)QQK",
                    "VNNIPS(ph)QSTR",
                    "SEGGGTGESS(ph)QGGTSK",
                    "GSLGS(ph)QGAKDEPEEELQK",
                    "GYGGS(ph)QGGGR",
                    "EEEGPAGEAAAS(ph)QPQAPTSVPGAR",
                    "EESGTIFGS(ph)QIK",
                    "VEVPSSAS(ph)QAK",
                    "SESPSLT(ph)QER",
                    "NVS(ph)QESLETK",
                ]
            )
        ),
        "Kinase Families": "ATR",
    }
)

CHECK1_seeds = pd.DataFrame(
    {
        "Modified sequence": sorted(
            set(
                [
                    "SLS(ph)GSPLK",
                    "FTS(ph)QAEK",
                    "KLS(ph)AEELER",
                    "SPS(ph)LLQSGAK",
                    "MLS(ph)FQGLAELAHR",
                    "SKS(ph)ATNLGR",
                    "AAS(ph)AIYR",
                    "RDS(ph)FDNCSLGESSK",
                ]
            )
        ),
        "Kinase Families": "sCHEK1",
    }
)

GSK3_seeds = pd.DataFrame(
    {
        "Modified sequence": sorted(
            set(
                [
                    "GGT(ph)PAGS(ph)ARGS(ph)PTRPNPPVR",
                    "GGT(ph)PAGS(ph)ARGS(ph)PTRPNPPVR",
                    "GGT(ph)PAGS(ph)ARGS(ph)PTRPNPPVR",
                    "TVT(ph)PASS(ph)AKTS(ph)PAK",
                    "TVT(ph)PASS(ph)AKTS(ph)PAK",
                    "TVT(ph)PASS(ph)AKTS(ph)PAK",
                    "SKS(ph)ATTT(ph)PSGS(ph)PR",
                    "SKS(ph)ATTT(ph)PSGS(ph)PR",
                    "SKS(ph)ATTT(ph)PSGS(ph)PR",
                    "SGGLQT(ph)PECLS(ph)REGS(ph)PIPHDPEFGSK",
                    "SGGLQT(ph)PECLS(ph)REGS(ph)PIPHDPEFGSK",
                    "SGGLQT(ph)PECLS(ph)REGS(ph)PIPHDPEFGSK",
                    "DLSTS(ph)PKPS(ph)PIPS(ph)PVLGR",
                    "DLSTS(ph)PKPS(ph)PIPS(ph)PVLGR",
                    "DLSTS(ph)PKPS(ph)PIPS(ph)PVLGR",
                    "AKTDHGAEIVYKS(ph)PVVS(ph)GDTS(ph)PR",
                    "AKTDHGAEIVYKS(ph)PVVS(ph)GDTS(ph)PR",
                    "AKTDHGAEIVYKS(ph)PVVS(ph)GDTS(ph)PR",
                    "PAS(ph)GPS(ph)SRPT(ph)SPCGK",
                    "PAS(ph)GPS(ph)SRPT(ph)SPCGK",
                    "PAS(ph)GPS(ph)SRPT(ph)SPCGK",
                    "EAS(ph)RESS(ph)RDTS(ph)PVR",
                    "EAS(ph)RESS(ph)RDTS(ph)PVR",
                    "EAS(ph)RESS(ph)RDTS(ph)PVR",
                    "IPRPSVS(ph)QGCS(ph)R",
                    "IPRPSVS(ph)QGCS(ph)R",
                    "MAS(ph)PPPS(ph)GPPS(ph)ATHT(ph)PFHQS(ph)PVEEK",
                    "MAS(ph)PPPS(ph)GPPS(ph)ATHT(ph)PFHQS(ph)PVEEK",
                    "MAS(ph)PPPS(ph)GPPS(ph)ATHT(ph)PFHQS(ph)PVEEK",
                    "MAS(ph)PPPS(ph)GPPS(ph)ATHT(ph)PFHQS(ph)PVEEK",
                    "KFELLPT(ph)PPLS(ph)PSR",
                    "S(ph)RPTS(ph)FADELAAR",
                    "HVT(ph)LPSS(ph)PR",
                    "STPQPPS(ph)GKTT(ph)PNSGDVQVTEDAVR",
                    "VLDTSSLTQS(ph)APAS(ph)PTNK",
                    "S(ph)RTPS(ph)ASNDDQQE",
                    "S(ph)RTPS(ph)ASNDDQQE",
                    "ADLVESLCSES(ph)TATS(ph)PV",
                    "DTYSDRS(ph)GSSS(ph)PDSEITELK",
                    "LAGGQTS(ph)QPTT(ph)PLTS(ph)PQR",
                    "LAGGQTS(ph)QPTT(ph)PLTS(ph)PQR",
                    "GT(ph)PSQS(ph)PVVGR",
                    "SQS(ph)AAVT(ph)PSSTTSSTR",
                    "KPSGDS(ph)QPSS(ph)PR",
                    "KPSGDS(ph)QPSS(ph)PR",
                    "VSPPAPGS(ph)APET(ph)PEDK",
                    "AASSS(ph)SPGS(ph)PVASS(ph)PSR",
                    "AASSS(ph)SPGS(ph)PVASS(ph)PSR",
                    "AASSS(ph)SPGS(ph)PVASS(ph)PSR",
                    "ELVGPPLAETVFT(ph)PKTS(ph)PENVQDR",
                    "AIS(ph)APTS(ph)PTR",
                    "AIS(ph)APTS(ph)PTR",
                    "GYYSPYSVS(ph)GSGS(ph)TAGSR",
                    "GS(ph)PATS(ph)PHLGR",
                    "TSPVVAPTS(ph)EPSS(ph)PLHTQLLK",
                    "S(ph)APSS(ph)PTLDCEK",
                    "GTEPS(ph)PGGT(ph)PQPSRPVS(ph)PAGPPEGVPEEAQPPR",
                    "GTEPS(ph)PGGT(ph)PQPSRPVS(ph)PAGPPEGVPEEAQPPR",
                    "GTEPS(ph)PGGT(ph)PQPSRPVS(ph)PAGPPEGVPEEAQPPR",
                    "LDQPVS(ph)APPS(ph)PR",
                    "LDQPVS(ph)APPS(ph)PR",
                    "HSVTAAT(ph)PPPS(ph)PTSGESGDLLSNLLQSPSSAK",
                    "MVSQS(ph)QPGS(ph)R",
                    "MVSQS(ph)QPGS(ph)R",
                    "STS(ph)TPTS(ph)PGPR",
                    "STS(ph)TPTS(ph)PGPR",
                    "VVSQS(ph)QPGS(ph)R",
                    "VVSQS(ph)QPGS(ph)R",
                    "LAES(ph)REQS(ph)PR",
                    "LAES(ph)REQS(ph)PR",
                    "SST(ph)PLHS(ph)PSPIR",
                    "EIPS(ph)ATQS(ph)PISK",
                    "GRLPNNSS(ph)RPST(ph)PTINVLESK",
                    "DGPSS(ph)APAT(ph)PTK",
                    "QFLIS(ph)PPAS(ph)PPVGWK",
                    "QFLIS(ph)PPAS(ph)PPVGWK",
                    "EKT(ph)PATT(ph)PEAR",
                    "SLLGDSAPTLHLNKGT(ph)PSQS(ph)PVVGR",
                    "SSESPSS(ph)SPSS(ph)PAR",
                    "SSESPSS(ph)SPSS(ph)PAR",
                    "AAAGPLDMS(ph)LPST(ph)PDIK",
                    "AAAGPLDMS(ph)LPST(ph)PDIK",
                    "YSQS(ph)APGS(ph)PVSAQPVIMAVPPRPSSLVAK",
                    "S(ph)APAS(ph)PTHPGLMS(ph)PR",
                    "S(ph)APAS(ph)PTHPGLMS(ph)PR",
                    "ETES(ph)APGS(ph)PR",
                    "ETES(ph)APGS(ph)PR",
                    "GPS(ph)QATS(ph)PIR",
                    "GPS(ph)QATS(ph)PIR",
                    "TPNNVVST(ph)PAPS(ph)PDASQLASSLSSQK",
                    "SIS(ph)REPS(ph)PALGPNLDGSGLLPR",
                    "METVSNASSSS(ph)NPSS(ph)PGR",
                    "METVSNASSSS(ph)NPSS(ph)PGR",
                    "SSVQGASS(ph)REGS(ph)PAR",
                    "SSVQGASS(ph)REGS(ph)PAR",
                    "APS(ph)RKDS(ph)LESDSSTAIIPHELIR",
                    "GYYSPYSVSGSGS(ph)TAGS(ph)RTGS(ph)R",
                    "GYYSPYSVSGSGS(ph)TAGS(ph)RTGS(ph)R",
                    "GYYSPYSVSGSGS(ph)TAGS(ph)RTGS(ph)R",
                    "TVDSQGPT(ph)PVCT(ph)PTFLER",
                    "TS(ph)RDTS(ph)PSSGSAVSSSK",
                    "ENS(ph)PAVS(ph)PTTNSTAPFGLKPR",
                    "RYPSS(ph)ISSS(ph)PQK",
                    "RYPSS(ph)ISSS(ph)PQK",
                    "SGAQASS(ph)TPLS(ph)PTR",
                    "SGAQASS(ph)TPLS(ph)PTR",
                    "ASVS(ph)GPNS(ph)PSETRR",
                ]
            )
        ),
        "Kinase Families": "GSK3",
    }
)

CDK9_seeds = pd.DataFrame(
    {
        "Modified sequence": sorted(
            set(
                [
                    "DVTNFTVGGFAPMS(ph)PR",
                    "TRT(ph)PAS(ph)INATPANINLADLTR",
                    "DHYQDPVPGIT(ph)PSSSSR",
                    "T(ph)PMYGSQTPLQDGSR",
                    "TPMYGSQT(ph)PLQDGSR",
                    "ETSSAT(ph)PGR",
                    "DPDANWDS(ph)PSR",
                    "SKS(ph)PRDPDANWDSPSR",
                    "YPESNRT(ph)PVKPSSVEEEDSFFR",
                    "VS(ph)PASSVDSNIPSSQGYK",
                    "MFS(ph)PMEEK",
                    "S(ph)PAQSDSTTQR",
                    "ELLS(ph)PLSEPDDR",
                    "QMSSQNS(ph)PSR",
                    "HRMS(ph)PGVAGS(ph)PR",
                    "HRMS(ph)PGVAGS(ph)PR",
                    "IST(ph)PQTNTVPIKPLIS(ph)TPPVSSQPK",
                    "IST(ph)PQTNTVPIKPLIS(ph)TPPVSSQPK",
                    "ENGFSS(ph)PPQIKDEPEDDGYFVPPK",
                    "NYGSPLISGST(ph)PK",
                    "ES(ph)PGAAATSSSGPQAQQHR",
                    "CAPSAGS(ph)PAAAVGR",
                    "YSPTS(ph)PTYSPTS(ph)PVYTPT(ph)SPK",
                    "YSPTS(ph)PTYSPTS(ph)PVYTPT(ph)SPK",
                    "YSPTS(ph)PTYSPTS(ph)PVYTPT(ph)SPK",
                    "SEPFS(ph)PSLRPEPPK",
                    "S(ph)PFEHSVEHK",
                    "NEEDEGHSNSS(ph)PRHS(ph)EAATAQR",
                    "HSS(ph)PHISR",
                ]
            )
        ),
        "Kinase Families": "sCDK9",
    }
)

CDK1213_seeds = pd.DataFrame(
    {
        "Modified sequence": sorted(
            set(
                [
                    "DTPGHGSGWAET(ph)PR",
                    "GDT(ph)PGHAT(ph)PGHGGATSSAR",
                    "GDT(ph)PGHAT(ph)PGHGGATSSAR",
                    "GGDSIGET(ph)PT(ph)PGASK",
                    "GGDSIGET(ph)PT(ph)PGASK",
                    "GSET(ph)PGATPGSK",
                    "GSETPGAT(ph)PGSK",
                    "IEEAMDGSET(ph)PQLFTVLPEK",
                    "IWDPTPSHTPAGAAT(ph)PGR",
                    "LSSWDQAET(ph)PGHT(ph)PSLR",
                    "LSSWDQAET(ph)PGHT(ph)PSLR",
                    "TMIIS(ph)PER",
                    "VLPPPAGYVPIRT(ph)PAR",
                    "WDET(ph)PGR",
                    "WDQTADQT(ph)PGAT(ph)PK",
                    "WDQTADQT(ph)PGAT(ph)PK",
                    "VVNGAAASQPPS(ph)KR",
                    "ILELT(ph)PEPDRPR",
                    "QTDPST(ph)PQQESSKPLGGIQPSSQTIQPK",
                    "TKPLT(ph)PSIGAK",
                    "T(ph)PTMPQEEAAEK",
                    "RT(ph)PTMPQEEAAACPPHILPPEK",
                    "HLLTDLPLPPELPGGDLS(ph)PPDS(ph)PEPK",
                    "HLLTDLPLPPELPGGDLS(ph)PPDS(ph)PEPK",
                    "ESKGS(ph)PVFLPR",
                    "SSGTAS(ph)SVAFT(ph)PLQGLEIVNPQAAEK",
                    "QETVADFT(ph)PK",
                ]
            )
        ),
        "Kinase Families": "sCDK1213",
    }
)

LIMK_seeds = pd.DataFrame(
    {
        "Modified sequence": sorted(set(["AS(ph)GVAVSDGVIK", "AS(ph)GVTVNDEVIK"])),
        "Kinase Families": "sLIMK",
    }
)

AAK1_seeds = pd.DataFrame(
    {
        "Modified sequence": sorted(
            set(
                [
                    "EQGGSGLGSGSS(ph)GGGGSTSGLGSGYIGR",
                    "SKS(ph)ATTTPSGSPR",
                    "ILSDVTHS(ph)AVFGVPASK",
                    "QGKGT(ph)ADETSK",
                    "VLFDNT(ph)GR",
                ]
            )
        ),
        "Kinase Families": "sAAK1",
    }
)

MAP2K47_seeds = pd.DataFrame(
    {
        "Modified sequence": sorted(
            set(
                [
                    "TAGTSFMMT(ph)PY(ph)VVTR",
                    "TACTNFMMT(ph)PY(ph)VVTR",
                    "TAGTSFMMTPY(ph)VVTR",
                    "TACTNFMMTPY(ph)VVTR",
                    "TAGTSFMMT(ph)PYVVTR",
                    "TACTNFMMT(ph)PYVVTR",
                ]
            )
        ),
        "Kinase Families": "sMAP2K47",
    }
)

MAP2K36_seeds = pd.DataFrame(
    {
        "Modified sequence": sorted(
            set(
                [
                    "HTDDEMT(ph)GY(ph)VATR",
                    "QADSEMT(ph)GY(ph)VVTR",
                    "HADAEMT(ph)GY(ph)VVTR",
                    "HTDDEMTGY(ph)VATR",
                    "QADSEMTGY(ph)VVTR",
                    "HADAEMTGY(ph)VVTR",
                    "HTDDEMT(ph)GYVATR",
                    "QADSEMT(ph)GYVVTR",
                    "HADAEMT(ph)GYVVTR",
                ]
            )
        ),
        "Kinase Families": "sMAP2K36",
    }
)

MAP2K12_seeds = pd.DataFrame(
    {
        "Modified sequence": sorted(
            set(
                [
                    "VADPDHDHTGFLT(ph)EY(ph)VATR",
                    "VADPDHDHTGFLTEY(ph)VATR",
                    "VADPDHDHTGFLT(ph)EYVATR",
                    "IADPEHDHTGFLT(ph)EY(ph)VATR",
                    "IADPEHDHTGFLT(ph)EYVATR",
                    "IADPEHDHTGFLTEY(ph)VATR",
                ]
            )
        ),
        "Kinase Families": "sMAP2K12",
    }
)

VALIDATED_KINASES = [
    "JNK",
    "MTOR",
    "ERK",
    "CDK2",
    "ATM",
    "CDK9/12/13",
    "sMAP2K36",
    "sMAP2K12",
    "PLK1",
    "p38",
    "MAPKAPKs",
    "sCHEK1",
    "RSKs",
    "AKTs",
    "MSKs",
    "TBK1",
    "CAMK2s",
    "ROCK1/2",
    "CDK1",
    "GSK3",
    "S6Ks",
    "AURKB",
    "CDK4/6",
    "sAAK1",
    "SRCs",
    "ATR",
    "CHEK2",
    "FAKs",
    # "RIPK2",
    # "FRK",
    # "HERs",
    # "SGKs",
    # "PAK1/2/3",
    # "MAP3K20",
    # "CDK12/13",
    # "CDK9",
    # "sATR",
    # "sGSK3",
    # "sCDK1213",
    # "sCDK9",
    # "EPHAs",
]

META_COLS = [
    "Patient_Identifier",
    "Program",
    "code_oncotree",
    "breadcrumb_oncotree",
    "tissue_topology",
]


def calculate_cytoplasmic_kinase_scores(
    results_folder: str,
    metadata_file: str,
    topas_kinase_substrate_file: str,
    expression_corrected_input: bool = False,
):
    results_folder = Path(results_folder)
    file_suffix = ""
    if expression_corrected_input:
        file_suffix = "_expressioncorrected"
    kinase_score_file = (
        results_folder
        / "topas_scores"
        / f"ck_substrate_phosphorylation_scores{file_suffix}.tsv"
    )
    if kinase_score_file.is_file():
        logger.info(
            f"Cytoplasmic kinase scoring skipped - found files already processed"
        )
        return

    logger.info("Construct joint modified sequence groups between cohort and decryptM")
    peptidoform_groups = get_joint_modified_sequence_groups(
        results_folder, topas_kinase_substrate_file
    )

    logger.info("Aggregate modified sequence groups in patient data")
    phospho = load_phospho_data(results_folder, file_suffix)
    pp_intensities_df = aggregate_patient_modified_sequence_groups(
        phospho, peptidoform_groups
    )

    logger.info("Aggregate modified sequence groups in decryptM data")
    automated_sites = load_confident_relationships(topas_kinase_substrate_file)
    automated_sites = aggregate_decryptm_modified_sequence_groups(
        automated_sites, peptidoform_groups
    )

    decryptM_kinases = get_patient_annotated_sites(pp_intensities_df, automated_sites)

    substrate_file = (
        results_folder / "topas_scores" / "ck_substrate_peptide_intensities.tsv"
    )
    write_substrate_peptides(
        pp_intensities_df, decryptM_kinases, substrate_file, kinases=VALIDATED_KINASES
    )

    logger.info("Compute cytoplasmic kinase scores")
    scores = compute_substrate_phosphorylation_scores(
        pp_intensities_df, decryptM_kinases, kinases=VALIDATED_KINASES
    )

    save_scores(scores, kinase_score_file)
    save_scores(scores, kinase_score_file, metadata_file)


def save_scores(
    scores: pd.DataFrame, kinase_score_file: Path, metadata_file: Optional[str] = None
):
    if metadata_file:
        metadata_df = sample_metadata.load(metadata_file)
        scores = merge_scores_with_sample_metadata(scores, metadata_df)
        kinase_score_file = kinase_score_file.with_name(
            kinase_score_file.stem + "_with_metadata.tsv"
        )

    scores.index.name = "Sample name"

    logger.info(f"Writing results to {kinase_score_file}")
    Path(kinase_score_file).parent.mkdir(exist_ok=True)

    scores.to_csv(kinase_score_file, sep="\t", float_format="%.4f")


def merge_scores_with_sample_metadata(
    scores: pd.DataFrame, metadata_df: pd.DataFrame
) -> pd.DataFrame:
    """Merge scores DataFrame with metadata on sample name."""
    score_cols = list(scores.columns)
    metadata_df["Sample name"] = "pat_" + metadata_df["Sample name"]
    merged = scores.merge(
        right=metadata_df.set_index("Sample name")[META_COLS],
        left_index=True,
        right_index=True,
        how="left",
    )
    return merged[META_COLS + score_cols]


def get_joint_modified_sequence_groups(
    results_folder: str, topas_kinase_substrate_file: str
):
    df_patients = read_cohort_modified_sequence_groups(results_folder)
    df_decryptM = get_decryptm_modified_sequence_groups(topas_kinase_substrate_file)
    df_annot = get_decryptm_modified_sequence_groups(
        topas_kinase_substrate_file, filter_for_confident_relationships=True
    )

    df_patients["Data set"] = "Patients"
    df_decryptM["Data set"] = "decryptM"
    df_annot["Data set"] = "Annotated"
    df_combined = pd.concat([df_patients, df_decryptM, df_annot])
    df_combined = explode_modified_sequence_groups(df_combined)
    df_combined

    # Aggregate duplicates
    df_combined = df_combined.groupby(["Modified sequence"])["Data set"].apply(set)
    df_combined = df_combined.reset_index()
    df_combined["dummy"] = 1  # dummy column for pa.aggregateModifiedSequenceGroups

    # Make the groups
    df_combined = pa.aggregateModifiedSequenceGroups(
        df_combined,
        experiment_cols=["dummy"],
        agg_cols={"Data set": lambda row: set.union(*row)},
        match_tolerance=2,
    )
    df_combined = df_combined[["Modified sequence group", "Data set"]]

    # Aggregate new Modified sequence group
    df_combined = df_combined.explode("Data set")
    df_combined.to_csv(
        f"{results_folder}/patient_decryptM_groups.txt", sep="\t", index=False
    )

    peptidoform_groups = df_combined[["Modified sequence group"]]
    peptidoform_groups["Modified sequence"] = peptidoform_groups[
        "Modified sequence group"
    ].str.split(";")
    peptidoform_groups = peptidoform_groups.explode("Modified sequence")
    return peptidoform_groups


def read_cohort_modified_sequence_groups(results_folder: str) -> pd.DataFrame:
    df_patients = pd.read_csv(
        results_folder / "preprocessed_pp2_agg_batchcorrected.csv",
        usecols=["Gene names", "Modified sequence group"],
    )
    return df_patients[["Modified sequence group"]]


def get_decryptm_modified_sequence_groups(
    topas_kinase_substrate_file: str, filter_for_confident_relationships: bool = False
) -> pd.DataFrame:
    df_decryptM = pd.read_csv(topas_kinase_substrate_file, sep="\t")
    df_decryptM = df_decryptM.rename(
        columns={
            "Modified sequence": "Modified sequence group",
            "Gene Names": "Gene names",
        }
    )
    if filter_for_confident_relationships:
        df_decryptM = df_decryptM[df_decryptM["Confident Relationship"]]
    df_decryptM = pd.concat(
        [df_decryptM, ATR_seeds, CHECK1_seeds, AAK1_seeds, GSK3_seeds]
    )
    df_decryptM = (
        df_decryptM.groupby(["Modified sequence group"])
        .first()
        .reset_index()[["Modified sequence group"]]
    )
    return df_decryptM


def explode_modified_sequence_groups(df: pd.DataFrame) -> pd.DataFrame:
    df["Modified sequence"] = df["Modified sequence group"].str.split(";")
    df = df.explode("Modified sequence")
    return df


def load_phospho_data(results_folder: str, file_suffix: str = "") -> pd.DataFrame:
    phospho = pd.read_csv(
        results_folder / f"preprocessed_pp2_agg_batchcorrected{file_suffix}.csv",
        index_col=[0, 1],
    )
    return phospho


def aggregate_patient_modified_sequence_groups(
    phospho: pd.DataFrame, peptidoform_groups: pd.DataFrame
) -> pd.DataFrame:
    phospho = 10**phospho
    phospho = phospho.reset_index()
    phospho["Modified sequence"] = phospho["Modified sequence group"].str.split(";")
    phospho = phospho.explode("Modified sequence")
    phospho = (
        phospho.drop(columns="Modified sequence group")
        .merge(right=peptidoform_groups, on="Modified sequence", how="left")
        .drop(columns="Modified sequence")
    )
    phospho = phospho.groupby(["Modified sequence group", "Gene names"]).mean()
    phospho = np.log10(phospho)
    return phospho


def load_confident_relationships(topas_kinase_substrate_file: str) -> pd.DataFrame:
    # Load annotations
    automated_sites = pd.read_csv(topas_kinase_substrate_file, sep="\t")
    automated_sites = automated_sites[automated_sites["Confident Relationship"]]
    automated_sites["Kinase Families"] = (
        automated_sites["Kinase Families"]
        .str.replace("CDK9", "CDK9/12/13")
        .str.replace("CDK12/13", "CDK9/12/13")
    )

    automated_sites = pd.concat(
        [
            automated_sites,
            ATR_seeds,
            CHECK1_seeds,
            GSK3_seeds,
            CDK9_seeds,
            CDK1213_seeds,
            AAK1_seeds,
            MAP2K47_seeds,
            MAP2K36_seeds,
            MAP2K12_seeds,
        ]
    )
    return automated_sites


def aggregate_decryptm_modified_sequence_groups(
    automated_sites: pd.DataFrame, peptidoform_groups: pd.DataFrame
):
    # Regroup based on all peptidoforms
    automated_sites["Modified sequence"] = automated_sites[
        "Modified sequence"
    ].str.split(";")
    automated_sites = automated_sites.explode("Modified sequence")
    automated_sites = automated_sites.merge(
        right=peptidoform_groups, on="Modified sequence", how="left"
    ).drop(columns="Modified sequence")
    automated_sites = (
        automated_sites.groupby("Modified sequence group")["Kinase Families"]
        .apply(set)
        .apply(sorted)
        .str.join(";")
    )
    return automated_sites


def get_patient_annotated_sites(
    pp_agg_df: pd.DataFrame, automated_sites: pd.DataFrame
) -> pd.Series:
    patient_modified_sequence_groups = pp_agg_df.reset_index()[
        pp_agg_df.index.names
    ].copy()
    decryptM_kinases = patient_modified_sequence_groups.merge(
        automated_sites, on="Modified sequence group", how="left"
    ).replace(np.nan, "")
    decryptM_kinases = decryptM_kinases.set_index(
        ["Modified sequence group", "Gene names"]
    )["Kinase Families"]
    return decryptM_kinases


def write_substrate_peptides(
    pp_intensities_df: pd.DataFrame,
    kinase_substrate_annotation_df: pd.Series,
    substrate_file: str,
    kinases: list[str] = None,
) -> pd.DataFrame:
    exploded_substrates_df = explode_series(kinase_substrate_annotation_df)
    if not kinases:
        kinases = exploded_substrates_df.unique()

    substrate_modified_sequence_groups = (
        exploded_substrates_df[exploded_substrates_df.isin(kinases)]
        .index.get_level_values("Modified sequence group")
        .drop_duplicates()
    )

    substrate_intensities_df = (
        kinase_substrate_annotation_df[substrate_modified_sequence_groups]
        .to_frame()
        .join(pp_intensities_df.loc[substrate_modified_sequence_groups], how="inner")
    )

    logger.info(f"Writing substrate intensities to {substrate_file}")
    Path(substrate_file).parent.mkdir(exist_ok=True)
    substrate_intensities_df.to_csv(substrate_file, sep="\t", float_format="%.4g")


def compute_substrate_phosphorylation_scores(
    pp_intensities_df: pd.DataFrame,
    kinase_substrate_annotation_df: pd.Series,
    kinases: list[str] = None,
    explode: bool = True,
) -> pd.DataFrame:
    if explode:
        kinase_substrate_annotation_df = explode_series(kinase_substrate_annotation_df)

    if kinases:
        kinase_substrate_annotation_df = kinase_substrate_annotation_df[
            kinase_substrate_annotation_df.isin(kinases)
        ]

    scores = pp_intensities_df.join(kinase_substrate_annotation_df, how="inner")
    scores = (
        scores.groupby(by=kinase_substrate_annotation_df.name)
        .progress_apply(
            z_aggregate,
            robust=False,
            standardize_input=False,
            standardize_output=True,
            agg_f="sum",
            clip_input=(-np.inf, np.inf),
            # clip_output=(-4, 4),
        )
        .T
    )
    return scores


def explode_series(s: pd.Series, delimiter: str = ";") -> pd.Series:
    s: pd.Series = s.str.split(delimiter)
    return s.explode()


def mad(x, sigma_scaling=1.4826017):
    if type(x) is pd.DataFrame:
        return (x.T - x.median(axis=1)).T.abs().median(axis=1) * sigma_scaling
    elif type(x) is pd.Series:
        return (x - x.median()).T.abs().median() * sigma_scaling


def z_aggregate(
    vals: pd.DataFrame,
    robust=False,
    standardize_input=True,
    center_output=True,
    standardize_output=True,
    agg_f="mean",
    clip_input=(-3, np.inf),
    clip_output=(-np.inf, np.inf),
) -> pd.Series:
    patient_vals = vals.filter(like="pat_")

    # Scale input (dataframe) in mean or robust way
    center, scale = get_center_and_scale(patient_vals, standardize_input, robust)
    z_vals = ((vals.T - center) / scale).T

    # Clip input to prevent destructive loss of a substrate for what ever reason
    z_vals = np.clip(z_vals, a_min=clip_input[0], a_max=clip_input[1])

    # Aggregate
    if agg_f == "mean":  # should be the best option
        agg_vals = z_vals.mean()
    elif agg_f == "median":
        agg_vals = z_vals.median()
    elif agg_f == "sum":  # sum is sensitive to unequal sizes (missing values)
        agg_vals = z_vals.sum()
    elif agg_f == "min":  # Min / Max is problematic in big aggregation sets
        agg_vals = z_vals.min()
    elif agg_f == "max":
        agg_vals = z_vals.max()
    else:
        raise ValueError("aggregation function is unknown...")

    agg_vals = agg_vals.replace(0, np.nan)

    # Scale output (series) in mean or robust
    if center_output or standardize_output:
        patient_agg_vals = agg_vals.filter(like="pat_")
        center, scale = get_center_and_scale(
            patient_agg_vals, standardize_output, robust, axis=0
        )
        agg_vals = (agg_vals - center) / scale

    # Clip output
    agg_vals = np.clip(agg_vals, a_min=clip_output[0], a_max=clip_output[1])

    return agg_vals


def get_center_and_scale(
    vals: pd.DataFrame, standardize: bool, robust: bool, axis: int = 1
):
    center = vals.median(axis=axis) if robust else vals.mean(axis=axis)
    scale = 1
    if standardize:
        scale = mad(vals) if robust else vals.std(axis=axis)
    return center, scale


"""
python3 -m topas_pipeline.topas.cytoplasmic_kinase_scoring -c config_patients.json [--expression-corrected]
"""
if __name__ == "__main__":
    import argparse

    from .. import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    parser.add_argument(
        "-e",
        "--expression-corrected",
        help="Use expression corrected phospho input.",
        action=argparse.BooleanOptionalAction,
    )
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    calculate_cytoplasmic_kinase_scores(
        results_folder=configs.results_folder,
        metadata_file=configs.metadata_annotation,
        topas_kinase_substrate_file=configs.clinic_proc.topas_kinase_substrate_file,
        expression_corrected_input=args.expression_corrected,
    )
