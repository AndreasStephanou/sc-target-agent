"""AnthropicModel is exercised against a fake SDK client so the parsing path is covered without a key."""
from types import SimpleNamespace

from sctarget.agents.model import AnthropicModel


class _Block(SimpleNamespace):
    def model_dump(self, exclude_none=True):
        return {k: v for k, v in vars(self).items() if v is not None}


class _FakeMessages:
    def __init__(self):
        self.last_kwargs = None

    def create(self, **kwargs):
        self.last_kwargs = kwargs
        return SimpleNamespace(
            content=[
                _Block(type="text", text="Loading first.", citations=None),
                _Block(type="tool_use", id="toolu_1", name="load_dataset", input={}),
            ],
            stop_reason="tool_use",
            usage=SimpleNamespace(input_tokens=120, output_tokens=15),
        )


def test_anthropic_model_parses_tool_use(monkeypatch):
    m = AnthropicModel(model="claude-test", api_key="sk-fake")
    fake = _FakeMessages()
    m.client = SimpleNamespace(messages=fake)
    resp = m.create(system="sys", messages=[{"role": "user", "content": "go"}], tools=[{"name": "load_dataset", "description": "x", "input_schema": {"type": "object", "properties": {}}}], role="analyst")
    assert resp.stop_reason == "tool_use"
    assert resp.text == "Loading first."
    assert len(resp.tool_uses) == 1 and resp.tool_uses[0].name == "load_dataset" and resp.tool_uses[0].id == "toolu_1"
    assert resp.usage == {"input_tokens": 120, "output_tokens": 15}
    assert resp.raw_content[1] == {"type": "tool_use", "id": "toolu_1", "name": "load_dataset", "input": {}}
    assert fake.last_kwargs["model"] == "claude-test" and fake.last_kwargs["tools"][0]["name"] == "load_dataset"
