from .frontend.data_structures import IdefixDataset, IdefixGrid, IdefixHierarchy
from .frontend.fields import IdefixFieldInfo
from .frontend.io import IdefixIOHandler

# register our dataset classes manually
from yt.utilities.object_registries import output_type_registry

output_type_registry["IdefixDataset"] = IdefixDataset
