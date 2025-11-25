# compute substrate phosphorylation scores (a.k.a. kinase activity) for cytoplasmic
# kinases using the confident substrates from Florian Bayer's decryptM experiments

import sys
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

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


def write_decryptM_annotations(
    topas_kinase_substrate_file: str, ck_substrate_annotations_file: str
):
    df_decryptM = pd.read_csv(
        topas_kinase_substrate_file,
        sep="\t",
        usecols=["Modified sequence", "Kinase Families", "Confident Relationship"],
    )
    df_decryptM = df_decryptM[df_decryptM["Confident Relationship"]]

    df_decryptM["Kinase Families"] = (
        df_decryptM["Kinase Families"]
        .str.replace("CDK9", "CDK9/12/13")
        .str.replace("CDK12/13", "CDK9/12/13")
    )
    df_decryptM = pd.concat(
        [
            df_decryptM,
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
    df_decryptM = df_decryptM.drop(columns="Confident Relationship")
    df_decryptM.to_csv(ck_substrate_annotations_file, sep="\t", index=False)


"""
python3 -m topas_pipeline.tools.create_ck_substrate_annotation_file -k <kinase_plausability_file> -c <ck_substrate_annotations_file>
"""
if __name__ == "__main__":
    import argparse

    from .. import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-k",
        "--kinase_plausability_file",
        required=True,
        help="Input path to Kinase plausability table from Flo's decryptM experiments.",
    )
    parser.add_argument(
        "-c",
        "--ck_substrate_annotations_file",
        required=True,
        help="Output path to ck substrate annotation file as input to the pipeline.",
    )
    args = parser.parse_args(sys.argv[1:])

    write_decryptM_annotations(
        topas_kinase_substrate_file=args.kinase_plausability_file,
        ck_substrate_annotations_file=args.ck_substrate_annotations_file,
    )
