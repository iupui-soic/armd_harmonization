"""Site loaders that build the canonical tables from staged Parquet."""

from .base import SiteLoader
from .dryad import ECUHLoader, StanfordLoader, UTSWLoader
from .mgb import MGBLoader

LOADERS: dict[str, type[SiteLoader]] = {
    "stanford": StanfordLoader,
    "mgb": MGBLoader,
    "utsw": UTSWLoader,
    "ecuh": ECUHLoader,
}


def get_loader(site: str, con=None) -> SiteLoader:
    try:
        return LOADERS[site](con)
    except KeyError:
        raise ValueError(f"no loader for site {site!r}") from None
