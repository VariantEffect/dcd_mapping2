import pytest
import pickle
from metadata_process import metadata_obtain
from blat_alignment import mave_to_blat
from transcript_selection import main

file = open("mappings.pickle", "rb")
mappings = pickle.load(file)
file.close()


@pytest.fixture(
    params=[
        "urn:mavedb:00000060-a-1",
        "urn:mavedb:00000041-a-1",
        "urn:mavedb:00000041-b-1",
    ]
)
def transcript_selection_dict(request):
    urn = request.param
    file_path = f"tests/data/{urn}"
    with open(file_path) as scoreset:
        mave_dat = metadata_obtain(scoreset)
    mave = mave_to_blat(mave_dat)
    tr = main(mave, mave_dat)
    return tr


def test_for_refseq_prot(transcript_selection_dict):
    mappings_dict = transcript_selection_dict
    computed_prot = mappings_dict["RefSeq_prot"]
    assert computed_prot == mappings[mappings_dict["urn"]][0]
    computed_nuc = mappings_dict["RefSeq_nuc"]
    assert computed_nuc == mappings[mappings_dict["urn"]][4]
