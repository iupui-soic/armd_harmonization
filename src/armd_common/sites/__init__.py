"""Site loaders that build the canonical tables from staged Parquet."""

from .base import SiteLoader
from .dryad import ECUHLoader, StanfordLoader, UTSWLoader
from .mgb import MGBLoader
from .mgb2 import MGB2Loader

LOADERS: dict[str, type[SiteLoader]] = {
    "stanford": StanfordLoader,
    "mgb": MGBLoader,
    "mgb2": MGB2Loader,
    "utsw": UTSWLoader,
    "ecuh": ECUHLoader,
}


def get_loader(site: str, con=None) -> SiteLoader:
    try:
        return LOADERS[site](con)
    except KeyError:
        raise ValueError(f"no loader for site {site!r}") from None
