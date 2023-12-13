import requests
import nest_asyncio
import asyncio
from cool_seq_tool.data_sources.mane_transcript_mappings import MANETranscriptMappings
from cool_seq_tool.data_sources.uta_database import UTADatabase
from Bio.Seq import Seq
from bs4 import BeautifulSoup
from mavedb_mapping.transcript_selection_helper import HelperFunctionsForBLATOutput
from mavedb_mapping import sr, dp

utadb = UTADatabase()
mane = MANETranscriptMappings()


nest_asyncio.apply()


async def mapq(locs: list, chrom: str):
    """
    Runs a query on UTADB to obtain transcripts
    Parameters
    ----------
        locs: list
           Locations list

        chrom: str
            Chromosome

        gsymb: str
            Gene symbol.


    Returns:
    --------
        list:
            List of transcripts from UTADB

    """
    transcript_lists = []
    for i in range(len(locs)):
        testquery = f"""select *
                            from uta_20210129.tx_exon_aln_v
                            where {locs[i][0]} between alt_start_i and alt_end_i
                            or {locs[i][1]} between alt_start_i and alt_end_i
                            and alt_ac = '{chrom}'"""

        out = await utadb.execute_query(testquery)
        tl = []
        for j in range(len(out)):
            if out[j]["tx_ac"].startswith("NR_") == False:
                tl.append(out[j]["tx_ac"])
        if tl != []:
            transcript_lists.append(tl)
    return transcript_lists

def using_uniprot(uniprot: str, ts: str):
    """
    Looks for transcripts using Uniprot ID

    Parameters
    ----------
        uniprot: str
            Uniprot ID

        ts: str
            Target Sequence


    Returns:
    --------
        start: str

        full_match:bool

    """
    try:
        url = "https://www.uniprot.org/uniprot/" + uniprot + ".xml"
        page = requests.get(url)
        page = BeautifulSoup(page.text)
        page = page.find_all("sequence")
        up = page[1].get_text()
        stri = str(ts)
        if up.find(stri) != -1:
            full_match = True
        else:
            full_match = False
        start = up.find(stri[:10])
        return start, full_match
    except:
        #no transcripts foundś
        return "NA", "NA"


def get_status(mane_trans: list):
    """
    Obtains status, RefSeq protein ID, and RefSeq nucleotide ID

    Parameters
    ----------
        mane_trans: list
            Transcripts obtained from MANE


    Returns:
    --------
        status: str
            status

        np: str
            RefSeq protein ID,

        nm: str
            RefSeq nucleotide ID


    """
    if len(mane_trans) == 1:
        np = mane_trans[0]["RefSeq_prot"]
        nm = mane_trans[0]["RefSeq_nuc"]
        status = "MANE Select"
    else:
        if mane_trans[0]["MANE_status"] == "MANE Select":
            np = mane_trans[0]["RefSeq_prot"]
            nm = mane_trans[0]["RefSeq_nuc"]
            status = "MANE Select"
        else:
            np = mane_trans[1]["RefSeq_prot"]
            nm = mane_trans[1]["RefSeq_nuc"]
            status = "MANE Plus Clinical"
    return status, np, nm


def from_mane_trans(dat: dict, mane_trans: list):
    """
    Obtains data for mapping using transcripts obtained from MANE

    Parameters
    ----------
        dat: dict
            Dictionary containing data from MaveDB scoresets

        mane_trans: list
            Transcripts obtained from MANE


    Returns:
    --------
        np: str
            RefSeq protein ID

        start: int

        full_match: bool

        nm: str
            RefSeq nucleotide ID

        status: str
            status

    """
    oseq = dat["target_sequence"]
    status, np, nm = get_status(mane_trans)
    if len(set(str(oseq))) > 4:
        stri = str(oseq)
    else:
        oseq = Seq(oseq)
        stri = str(oseq.translate(table=1)).replace("*", "")

    if str(sr[np]).find(stri) != -1:
        full_match = True
    else:
        full_match = False
    start = str(sr[np]).find(stri[:10])
    mappings_list = [np, start, dat["urn"], full_match, nm, status]
    return np, start, full_match, nm, status


async def np(nm):
    testquery = (
        f"SELECT pro_ac FROM uta_20210129.associated_accessions WHERE tx_ac = '{nm}'"
    )
    out = await utadb.execute_query(testquery)
    try:
        return out[0]["pro_ac"]
    except:
        return out


def no_mane_trans(isect, dat):
    """
    Obtains data for mapping if transcripts were not found using MANE

    Parameters
    ----------
        isect: list
            List with Refseq nucleotide identifiers

        dat: dict
            Dictionary containing data from MaveDB scoresets


    Returns:
    --------
        np: str
            RefSeq protein ID

        start: int

        full_match: bool

        nm: str
            RefSeq nucleotide ID

        status: str
            status

    """
    trans_lens = []
    for i in range(len(isect)):
        trans_lens.append(len(str(sr[isect[i]])))
        loc = trans_lens.index(max(trans_lens))
        nm = isect[loc]
        np = asyncio.run(np(nm))

    if np != []:
        oseq = dat["target_sequence"]

        if len(set(str(oseq))) > 4:
            stri = str(oseq)
        else:
            oseq = Seq(oseq)
            stri = str(oseq.translate(table=1)).replace("*", "")

        if str(sr[np]).find(stri) != -1:
            full_match = True
        else:
            full_match = False
        start = str(sr[np]).find(stri[:10])
        status = "Longest Compatible"
        return np, start, full_match, nm, status


def main(mave_blat_dict: dict, dat: dict) -> dict:
    """
    Returns dictionary after trancsript selection

    Parameters
    ----------

        mave_blat_dict: dict
            Dictionary obtained after BLAT Alignment

        dat: dict
            Dictionary containing data from MaveDB scoresets

    Returns
    -------
        mappings_dict:
            Dictionary after transcript selections"""

    helper = HelperFunctionsForBLATOutput(mave_blat_dict)

    if dat["target_type"] == "Protein coding" or dat["target_type"] == "protein_coding":
      
        if not helper.is_human():
            raise ValueError("Invalid Scoreset")

        locs = helper.get_locs_list()
        chrom = helper.get_chr(dp)
        ts = asyncio.run(mapq(locs, chrom))

        try:
            for i in ts:
                i = set(i)
            intersection_set = set.intersection(set(i))
            isect = list(intersection_set)

        except:
            start, full_match = using_uniprot(dat["uniprot_id"], dat["target_sequence"])
            np = nm = status = "NA"
            mappings_dict = {
                "urn": dat["urn"],
                "uniprot_id": dat["uniprot_id"],
                "start": start,
                "full_match": full_match,
                "RefSeq_prot": "NA",
                "RefSeq_nuc": "NA",
                "status": "NA",
                "gsymb": "NA",
            }
            return mappings_dict

        mane_trans = mane.get_mane_from_transcripts(isect)
        
        
        if mane_trans != []:
            #selecting correct transcript using chrom from transcripts obtained
            for i in mane_trans:
                if i['GRCh38_chr'] == chrom:
                    mane_trans = [i]
            gsymb = mane_trans[0]["symbol"]
            
            np, start, full_match, nm, status = from_mane_trans(dat, mane_trans)
        else:
            np, start, full_match, nm, status = no_mane_trans(isect, dat)

        mappings_dict = {
            "urn": dat["urn"],
            "uniprot_id": dat["uniprot_id"],
            "start": start,
            "full_match": full_match,
            "RefSeq_prot": np,
            "RefSeq_nuc": nm,
            "status": status,
            "gsymb": gsymb,
            "uniprot_id": dat["uniprot_id"],
        }
        return mappings_dict