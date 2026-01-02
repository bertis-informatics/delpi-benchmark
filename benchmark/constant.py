from pathlib import Path

from benchmark.tools.msgf import MSGFReader
from benchmark.tools.sage import SageReader
from benchmark.tools.diann import DIANNReader
from benchmark.tools.alphadia import AlphaDIAReader
from benchmark.tools.diabert import DIABertReader
from benchmark.tools.delpi import DelPiReader


tool_to_reader_class = {
    "msgf": MSGFReader,
    "sage": SageReader,
    "diann-1.8.1": DIANNReader,
    "alphadia": AlphaDIAReader,
    "diabert": DIABertReader,
    "delpi": DelPiReader,
}

tool_color_map = {
    "delpi": "#004E98",
    "diann-1.8.1": "#E94F37",
    "alphadia": "#F79256",
    "diabert": "#F9C74F",
    "sage": "#008B8B",
    "msgf": "#8CC64D",
}

tool_display_name_map = {
    "delpi": "DelPi",
    "diann-1.8.1": "DIA-NN",
    "alphadia": "AlphaDIA",
    "diabert": "DIA-BERT",
    "sage": "Sage",
    "msgf": "MS-GF+",
}
