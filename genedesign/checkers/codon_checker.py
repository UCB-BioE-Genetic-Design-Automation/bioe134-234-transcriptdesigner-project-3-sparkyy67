import sys
import csv
from collections import Counter

class CodonChecker:
    """
    Description:
    This class checks codon usage for a given CDS and calculates metrics like:
    1. Codon Diversity: Fraction of unique codons.
    2. Rare Codon Count: Number of rare codons in the sequence.
    3. Codon Adaptation Index (CAI): Based on codon usage frequencies.

    Input (run method):
    cds (List[str]): A list of codons representing the coding sequence (e.g., ['ATG', 'TAA', 'CGT']).

    Output:
    Tuple[bool, float, int, float]: A tuple containing:
        - codons_above_board (bool): True if the CDS passes the thresholds, False otherwise.
        - codon_diversity (float): Fraction of unique codons.
        - rare_codon_count (int): Number of rare codons.
        - cai_value (float): CAI value.
    """

    codon_frequencies: dict[str, float]
    rare_codons: list[str]
    rare_codon_threshold: float

    def initiate(self) -> None:
        """
        Loads codon usage data from a file and sets up the codon frequencies and rare codons.
        """
        codon_usage_file = 'genedesign/data/codon_usage.txt'
        self.codon_frequencies = {}
        self.rare_codons = []
        self.rare_codon_threshold = 0.1

        with open(codon_usage_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 4:
                    continue
                codon = row[0].strip()
                usage_freq = float(row[2].strip())
                self.codon_frequencies[codon] = usage_freq
                if usage_freq < self.rare_codon_threshold:
                    self.rare_codons.append(codon)

    def run(self, cds: list[str]) -> tuple[bool, float, int, float]:
        """
        Calculates codon diversity, rare codon count, and CAI for the provided CDS.

        :param cds: List of codons representing the CDS.
        :return: Tuple containing a boolean, codon diversity, rare codon count, and CAI score.
        """
        if not cds:
            return False, 0.0, 0, 0.0

        # Diversity = unique codons / total codons
        codon_counts = Counter(cds)
        codon_diversity = len(codon_counts) / len(cds)

        # Count rare codons
        rare_codon_count = sum(codon_counts[codon] for codon in self.rare_codons if codon in cds)

        # CAI as geometric mean of codon frequencies
        cai_numerators = [self.codon_frequencies.get(codon, 0.01) for codon in cds]
        cai_product = 1
        for freq in cai_numerators:
            cai_product *= freq
        cai_value = cai_product ** (1 / len(cai_numerators)) if cai_numerators else 0.0

        # Thresholds
        diversity_threshold = 0.3
        rare_codon_limit = 5
        cai_threshold = 0.15

        codons_above_board = (codon_diversity >= diversity_threshold and
                              rare_codon_count <= rare_codon_limit and
                              cai_value >= cai_threshold)

        return codons_above_board, codon_diversity, rare_codon_count, cai_value


if __name__ == "__main__":
    codon_checker = CodonChecker()
    codon_checker.initiate()
    cds = ['ATG', 'CAA', 'GGG', 'TAA']
    codons_above_board, codon_diversity, rare_codon_count, cai_value = codon_checker.run(cds)
    print(f"CDS: {cds}")
    print(f"Codons Above Board: {codons_above_board}")
    print(f"Codon Diversity: {codon_diversity}")
    print(f"Rare Codon Count: {rare_codon_count}")
    print(f"Codon Adaptation Index (CAI): {cai_value}")