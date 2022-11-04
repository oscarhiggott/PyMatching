import pathlib
import os

THIS_DIR = pathlib.Path(__file__).parent.resolve()
DATA_DIR = os.path.join(pathlib.Path(THIS_DIR).parent.parent.absolute(), "data")
