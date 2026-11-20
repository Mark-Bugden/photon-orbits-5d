"""Checks that tests/fixtures/orbits.json, which the JS port's tests replay, still matches the reference model."""

import json
import math

from photon_orbits import fixtures

# The fixtures are written on Windows and CI runs Linux, whose sin, cos and atan can differ in the last bit,
# so numbers are compared to within this relative tolerance. A change to the model shows up far above it.
TOLERANCE = 1e-10


def mismatches(committed, current, path="fixture"):
    """Every place where the two JSON structures differ, beyond TOLERANCE for numbers."""
    if isinstance(committed, dict) and isinstance(current, dict):
        if committed.keys() != current.keys():
            return [f"{path}: keys {sorted(committed)} vs {sorted(current)}"]
        return [m for key in committed for m in mismatches(committed[key], current[key], f"{path}.{key}")]
    if isinstance(committed, list) and isinstance(current, list):
        if len(committed) != len(current):
            return [f"{path}: length {len(committed)} vs {len(current)}"]
        return [m for i, (a, b) in enumerate(zip(committed, current)) for m in mismatches(a, b, f"{path}[{i}]")]
    if isinstance(committed, (int, float)) and isinstance(current, (int, float)) and not isinstance(committed, bool):
        return [] if math.isclose(committed, current, rel_tol=TOLERANCE, abs_tol=TOLERANCE) else [f"{path}: {committed} vs {current}"]
    return [] if committed == current else [f"{path}: {committed!r} vs {current!r}"]


def test_the_js_fixtures_are_current():
    committed = json.loads(fixtures.OUT.read_text(encoding="utf-8"))
    current = json.loads(json.dumps(fixtures.build()))
    found = mismatches(committed, current)
    assert not found, (
        "tests/fixtures/orbits.json is out of date: run `uv run python -m photon_orbits.fixtures`\n"
        + "\n".join(found[:10])
    )


def test_a_model_change_would_be_caught():
    current = json.loads(json.dumps(fixtures.build()))
    changed = json.loads(json.dumps(current))
    changed["orbits"][7]["table"]["phi"][5] *= 1 + 1e-8
    assert mismatches(changed, current)


def test_the_cases_cover_what_they_say():
    built = fixtures.build()
    assert all(case["reason"] for case in built["noOrbits"])
    assert {case["reason"] for case in built["noOrbits"]} == {"no band", "at horizon", "no horizon"}
    assert all(len(case["edges"]) == 96 for case in built["regions"])
    assert any(case["orbit"]["Phi"] == 0 and case["Phi"] != 0 for case in built["orbits"])       # a snapped momentum
