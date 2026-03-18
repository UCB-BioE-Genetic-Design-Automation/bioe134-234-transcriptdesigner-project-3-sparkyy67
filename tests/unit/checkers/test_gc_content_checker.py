import pytest
from genedesign.checkers.gc_content_checker import GCContentChecker

@pytest.fixture
def checker():
    c = GCContentChecker()
    c.initiate()
    return c

def test_normal_gc(checker):
    # 50% GC - should pass
    result, gc = checker.run("ATGCATGC")
    assert result == True
    assert abs(gc - 0.5) < 0.01

def test_too_low_gc(checker):
    # All AT - should fail
    result, gc = checker.run("AAAATTTT")
    assert result == False
    assert gc == 0.0

def test_too_high_gc(checker):
    # All GC - should fail
    result, gc = checker.run("GGGGCCCC")
    assert result == False
    assert gc == 1.0

def test_boundary_low(checker):
    # Exactly 40% GC - should pass
    dna = "GC" * 4 + "AT" * 6  # 8 GC out of 20 = 40%
    result, gc = checker.run(dna)
    assert result == True
    assert abs(gc - 0.4) < 0.01

def test_empty_sequence(checker):
    result, gc = checker.run("")
    assert result == False