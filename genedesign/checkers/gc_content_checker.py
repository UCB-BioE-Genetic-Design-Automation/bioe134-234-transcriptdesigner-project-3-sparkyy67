class GCContentChecker:
    """
    Checks that a DNA sequence has GC content between 40% and 60%.
    Sequences outside this range tend to have secondary structure or
    expression problems in E. coli.

    Input:
        dna (str): A DNA sequence string.

    Output:
        Tuple[bool, float]:
            - True if GC content is within bounds, False otherwise.
            - The GC content as a float (e.g. 0.52 for 52%).
    """

    def initiate(self) -> None:
        self.min_gc = 0.40
        self.max_gc = 0.60

    def run(self, dna: str) -> tuple[bool, float]:
        if not dna:
            return False, 0.0
        dna = dna.upper()
        gc_count = dna.count('G') + dna.count('C')
        gc_content = gc_count / len(dna)
        passed = self.min_gc <= gc_content <= self.max_gc
        return passed, gc_content


if __name__ == "__main__":
    checker = GCContentChecker()
    checker.initiate()
    result, gc = checker.run("ATGCGATCGATCG")
    print(f"Passed: {result}, GC Content: {gc:.2%}")