from .data_structures import IdefixDumpDataset, IdefixGrid, IdefixHierarchy
from .fields import IdefixFieldInfo
from .io import IdefixIOHandler

# register our dataset classes manually
from yt.utilities.object_registries import output_type_registry

output_type_registry["IdefixDumpDataset"] = IdefixDumpDataset
