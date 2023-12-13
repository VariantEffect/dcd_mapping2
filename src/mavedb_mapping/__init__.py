
import os
from biocommons.seqrepo import SeqRepo

from cool_seq_tool.data_sources.uta_database import UTADatabase

from ga4gh.vrs.dataproxy import SeqRepoDataProxy

from gene.query import QueryHandler
from gene.database import create_db

__version__ = "0.1.0"

path_to_seqrepo = os.getenv("PATH_TO_SEQREPO", "/usr/local/share/seqrepo/latest")
path_to_hg38_file = os.getenv("HG38_FILE", "hg38.2bit")  # default- in current directory
data_file_path = "tests/data/"

# utadb = UTADatabase(db_pwd="uta")
qh = QueryHandler(create_db())
sr = SeqRepo(path_to_seqrepo, writeable=True)
dp = SeqRepoDataProxy(sr=sr)