"""Tool registry: every function the agent may call is registered here with a JSON schema.

The agent never executes arbitrary analysis code against the data. It chooses from
these tools and sets their parameters. That is what makes a run auditable and
reproducible: the run log is a list of (tool, params, result) triples.
"""
from __future__ import annotations

import inspect
from dataclasses import dataclass, field
from typing import Any, Callable

ToolFn = Callable[..., dict[str, Any]]


@dataclass
class ToolSpec:
    name: str
    description: str
    input_schema: dict[str, Any]
    fn: ToolFn
    roles: set[str] = field(default_factory=lambda: {"*"})
    mutates_state: bool = False

    def anthropic_schema(self) -> dict[str, Any]:
        return {"name": self.name, "description": self.description, "input_schema": self.input_schema}

    def allowed_for(self, role: str) -> bool:
        return "*" in self.roles or role in self.roles


class ToolRegistry:
    def __init__(self) -> None:
        self._tools: dict[str, ToolSpec] = {}

    def register(self, spec: ToolSpec) -> None:
        if spec.name in self._tools:
            raise ValueError(f"tool {spec.name} already registered")
        self._tools[spec.name] = spec

    def get(self, name: str) -> ToolSpec:
        if name not in self._tools:
            raise KeyError(f"unknown tool: {name}")
        return self._tools[name]

    def names(self) -> list[str]:
        return sorted(self._tools)

    def for_role(self, role: str) -> list[ToolSpec]:
        return [t for t in self._tools.values() if t.allowed_for(role)]

    def schemas_for_role(self, role: str) -> list[dict[str, Any]]:
        return [t.anthropic_schema() for t in self.for_role(role)]

    def __contains__(self, name: str) -> bool:
        return name in self._tools

    def __len__(self) -> int:
        return len(self._tools)


REGISTRY = ToolRegistry()


def tool(
    name: str,
    description: str,
    input_schema: dict[str, Any],
    roles: set[str] | None = None,
    mutates_state: bool = False,
) -> Callable[[ToolFn], ToolFn]:
    """Decorator registering a function as an agent-callable tool.

    The decorated function must accept ``session`` as its first argument plus keyword
    arguments matching ``input_schema['properties']``.
    """

    def deco(fn: ToolFn) -> ToolFn:
        sig = inspect.signature(fn)
        params = list(sig.parameters)
        if not params or params[0] != "session":
            raise TypeError(f"tool {name}: first parameter must be 'session'")
        schema_keys = set(input_schema.get("properties", {}))
        fn_keys = set(params[1:])
        missing = schema_keys - fn_keys
        if missing:
            raise TypeError(f"tool {name}: schema properties {missing} not in function signature")
        input_schema.setdefault("type", "object")
        input_schema.setdefault("additionalProperties", False)
        REGISTRY.register(
            ToolSpec(
                name=name,
                description=description.strip(),
                input_schema=input_schema,
                fn=fn,
                roles=roles or {"*"},
                mutates_state=mutates_state,
            )
        )
        return fn

    return deco


def validate_args(spec: ToolSpec, args: dict[str, Any]) -> dict[str, Any]:
    """Light-weight validation: required keys present, no unknown keys, basic types."""
    props = spec.input_schema.get("properties", {})
    required = spec.input_schema.get("required", [])
    unknown = set(args) - set(props)
    if unknown:
        raise ValueError(f"{spec.name}: unknown arguments {sorted(unknown)}")
    for r in required:
        if r not in args:
            raise ValueError(f"{spec.name}: missing required argument '{r}'")
    type_map = {"string": str, "integer": int, "number": (int, float), "boolean": bool, "array": list, "object": dict}
    for k, v in args.items():
        t = props[k].get("type")
        if t is None or v is None:
            continue
        if isinstance(t, list):
            ok = any(isinstance(v, type_map[tt]) for tt in t if tt in type_map) or ("null" in t and v is None)
        else:
            ok = isinstance(v, type_map.get(t, object))
        if t == "integer" and isinstance(v, bool):
            ok = False
        if not ok:
            raise ValueError(f"{spec.name}: argument '{k}' should be {t}, got {type(v).__name__}")
    return args
