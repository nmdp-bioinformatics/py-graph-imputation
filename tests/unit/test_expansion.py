# -*- coding: utf-8 -*-

#
#    grim Graph Imputation
#    Copyright (c) 2021 Be The Match operated by National Marrow Donor Program. All Rights Reserved.
#
#    This library is free software; you can redistribute it and/or modify it
#    under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation; either version 3 of the License, or (at
#    your option) any later version.
#
#    This library is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#    License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library;  if not, write to the Free Software Foundation,
#    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.
#
#    > http://www.fsf.org/licensing/licenses/lgpl.html
#    > http://www.opensource.org/licenses/lgpl-license.php
#

"""Tests for ambiguity (XX / MAC / serology) expansion in imputation_service.

These cover the py-ard integration layer that turns ambiguous GL strings into
'/'-separated allele lists the graph engine can impute over. The pure
detection/attachment helpers run without building a graph or touching py-ard;
the actual expansion test is skipped if py-ard (or its reference DB) is not
available in the environment.
"""

import pytest

import imputation_service as svc


# ---------------------------------------------------------------------------
# _has_ambiguous_codes: detects XX codes, MAC codes, and serology
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "gl_string,expected",
    [
        # Plain, unambiguous 2-field genotypes -> no expansion needed.
        ("A*01:01+A*02:01^B*07:02+B*08:01^DRB1*03:01+DRB1*15:01", False),
        # Explicit '/'-ambiguity is understood by the engine directly.
        ("A*01:01/A*01:02+A*02:01", False),
        ("", False),
        # XX codes.
        ("A*01:XX+A*02:01", True),
        ("A*01:xx+A*02:01", True),  # case-insensitive
        ("DRB1*15:XX+DRB1*03:01", True),
        # MAC codes (uppercase alpha second field).
        ("A*01:AB+A*02:01", True),
        ("B*07:BXYZ+B*08:01", True),
        # Serology.
        ("A1+A2", True),
        ("Cw7+B8", True),
    ],
)
def test_has_ambiguous_codes(gl_string, expected):
    assert svc._has_ambiguous_codes(gl_string) is expected


def test_has_ambiguous_codes_none():
    assert svc._has_ambiguous_codes(None) is False


# ---------------------------------------------------------------------------
# _attach_expansion: records how the input was expanded, on the response dict
# ---------------------------------------------------------------------------


def test_attach_expansion_when_expanded():
    resp = {"gl_string": "A*01:01/A*01:02+A*02:01"}
    svc.ImputationService._attach_expansion(
        resp,
        original_gl_string="A*01:XX+A*02:01",
        gl_expanded=True,
        expand_note=None,
    )
    assert resp["expanded"] is True
    assert resp["input_gl_string"] == "A*01:XX+A*02:01"
    assert "expand_note" not in resp


def test_attach_expansion_when_not_expanded():
    resp = {"gl_string": "A*01:01+A*02:01"}
    svc.ImputationService._attach_expansion(
        resp,
        original_gl_string="A*01:01+A*02:01",
        gl_expanded=False,
        expand_note=None,
    )
    # A pass-through carries no expansion metadata.
    assert "expanded" not in resp
    assert "input_gl_string" not in resp
    assert "expand_note" not in resp


def test_attach_expansion_records_note():
    resp = {"gl_string": "A*01:XX+A*02:01"}
    svc.ImputationService._attach_expansion(
        resp,
        original_gl_string="A*01:XX+A*02:01",
        gl_expanded=False,
        expand_note="py-ard not available; codes not expanded",
    )
    assert resp["expand_note"] == "py-ard not available; codes not expanded"
    assert "expanded" not in resp


# ---------------------------------------------------------------------------
# expand_gl_string: graceful fallback + real py-ard expansion
# ---------------------------------------------------------------------------


def test_expand_gl_string_falls_back_without_pyard(monkeypatch):
    """When py-ard is unavailable, the input is returned unchanged with a note."""
    monkeypatch.setattr(svc, "_get_ard", lambda imgt_version=None: None)
    gl = "A*01:XX+A*02:01"
    out, note = svc.expand_gl_string(gl)
    assert out == gl
    assert note is not None
    assert "py-ard" in note


def _pyard_available():
    ard = svc._get_ard()
    return ard is not None


@pytest.mark.skipif(
    not _pyard_available(),
    reason="py-ard (or its reference DB) not available in this environment",
)
def test_expand_gl_string_expands_xx_code():
    """A real XX code expands to a '/'-separated 2-field allele list."""
    out, note = svc.expand_gl_string("A*01:XX+A*02:01")
    assert note is None
    # XX expands into multiple members joined by '/', at 2-field resolution.
    assert "/" in out
    assert "A*01:01" in out
    # 2-field lgx reduction: no 3rd field like ':01:01' in the expansion.
    assert "A*01:XX" not in out
