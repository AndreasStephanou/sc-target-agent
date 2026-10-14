"""Importing this package registers every tool in :data:`sctarget.tools.registry.REGISTRY`."""
from . import analysis, knowledge, targets  # noqa: F401  (registration side effects)
from .registry import REGISTRY, ToolSpec, validate_args

__all__ = ["REGISTRY", "ToolSpec", "validate_args", "analysis", "knowledge", "targets"]
