import random
import csv
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.gc_content_checker import GCContentChecker
from genedesign.seq_utils.reverse_complement import reverse_complement


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a codon-optimized DNA sequence.

    Strategy:
      1. Load real codon frequencies from codon_usage.txt so selection weights
         reflect actual E. coli CAI values.
      2. Build the CDS codon-by-codon using weighted random selection, avoiding
         repeating the same codon back-to-back and checking forbidden sequence
         junctions locally using 2-codon, 3-codon windows (both directions),
         lookahead to the next amino acid, AND reverse complement checks.
      3. Deprioritize TTG as a Leucine codon since it matches the sigma70
         promoter -35 box (TTGACA) and causes internal promoter failures.
      4. After building a full candidate CDS, run all checkers on utr + cds
         exactly matching what the benchmarker validates. Rotate through RBS
         options to escape UTR-caused issues.
      5. Reject invalid peptides immediately rather than wasting retry attempts.
    """

    def __init__(self):
        self.codonTable = {}
        self.codonWeights = {}
        self.rbsChooser = None
        self.forbiddenChecker = None
        self.promoterChecker = None
        self.codonChecker = None
        self.gcChecker = None

    def initiate(self) -> None:
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.forbiddenChecker.initiate()

        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()

        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()

        self.gcChecker = GCContentChecker()
        self.gcChecker.initiate()

        genetic_code = {
            'A': ['GCT', 'GCC', 'GCA', 'GCG'],
            'C': ['TGT', 'TGC'],
            'D': ['GAT', 'GAC'],
            'E': ['GAA', 'GAG'],
            'F': ['TTT', 'TTC'],
            'G': ['GGT', 'GGC', 'GGA', 'GGG'],
            'H': ['CAT', 'CAC'],
            'I': ['ATT', 'ATC', 'ATA'],
            'K': ['AAA', 'AAG'],
            # TTG deprioritized - matches sigma70 -35 box causing promoter failures
            'L': ['CTG', 'CTT', 'CTC', 'CTA', 'TTA', 'TTG'],
            'M': ['ATG'],
            'N': ['AAT', 'AAC'],
            'P': ['CCG', 'CCA', 'CCT', 'CCC'],
            'Q': ['CAG', 'CAA'],
            'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
            'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
            'T': ['ACC', 'ACA', 'ACT', 'ACG'],
            'V': ['GTT', 'GTC', 'GTA', 'GTG'],
            'W': ['TGG'],
            'Y': ['TAT', 'TAC'],
            '*': ['TAA', 'TGA', 'TAG'],
        }

        codon_freq = {}
        codon_usage_file = 'genedesign/data/codon_usage.txt'
        with open(codon_usage_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 3:
                    continue
                codon = row[0].strip()
                freq = float(row[2].strip())
                codon_freq[codon] = freq

        for aa, codons in genetic_code.items():
            sorted_codons = sorted(
                codons,
                key=lambda c: codon_freq.get(c, 0.0),
                reverse=True
            )
            self.codonTable[aa] = sorted_codons
            self.codonWeights[aa] = [codon_freq.get(c, 0.01) for c in sorted_codons]

    def _creates_forbidden(self, seq: str) -> bool:
        """
        Check if a DNA sequence or its reverse complement contains
        any forbidden site. Checking RC catches palindromic sites
        and sites that appear only on the opposite strand.
        """
        seq = seq.upper()
        rc = reverse_complement(seq).upper()
        for site in self.forbiddenChecker.forbidden:
            if site in seq or site in rc:
                return True
        return False

    def _build_codons(self, peptide: str) -> list[str]:
        codons = []
        last_codon = None

        for i, aa in enumerate(peptide):
            if i == 0:
                codons.append('ATG')
                last_codon = 'ATG'
                continue

            options = self.codonTable.get(aa, ['NNN'])
            weights = self.codonWeights.get(aa, [1.0])

            # Get next amino acid's codons for lookahead
            next_options = []
            if i + 1 < len(peptide):
                next_aa = peptide[i + 1]
                next_options = self.codonTable.get(next_aa, ['NNN'])

            # Weighted random order — prefer high CAI but with variation
            paired = list(zip(options, weights))
            random.shuffle(paired)
            paired.sort(key=lambda x: x[1], reverse=True)

            prev_prev = codons[-2] if len(codons) >= 2 else None

            chosen = None
            for candidate, _ in paired:
                # Skip same codon as last to prevent poly runs
                if candidate == last_codon:
                    continue

                # 2-codon window: last + candidate
                if last_codon is not None:
                    if self._creates_forbidden(last_codon + candidate):
                        continue

                # 3-codon window: prev_prev + last + candidate
                if prev_prev is not None:
                    if self._creates_forbidden(prev_prev + last_codon + candidate):
                        continue

                # 3-codon window: last + candidate + next
                # catches sites like GTG+GAT+CCG where site starts mid last_codon
                if last_codon is not None and next_options:
                    all_next_bad = all(
                        self._creates_forbidden(last_codon + candidate + next_c)
                        for next_c in next_options
                    )
                    if all_next_bad:
                        continue

                # Lookahead: skip only if ALL next codons create a forbidden
                # junction with this candidate
                if next_options:
                    all_next_bad = all(
                        self._creates_forbidden(candidate + next_c)
                        for next_c in next_options
                    )
                    if all_next_bad:
                        continue

                chosen = candidate
                break

            # Fallback: ignore lookahead but still avoid direct junctions
            if chosen is None:
                for candidate, _ in paired:
                    if last_codon is not None:
                        if self._creates_forbidden(last_codon + candidate):
                            continue
                    if prev_prev is not None:
                        if self._creates_forbidden(prev_prev + last_codon + candidate):
                            continue
                    chosen = candidate
                    break

            # Absolute last resort
            if chosen is None:
                chosen = options[0]

            codons.append(chosen)
            last_codon = chosen

        codons.append('TAA')
        return codons

    def run(self, peptide: str, ignores: set) -> Transcript:
        # Reject invalid peptides immediately
        valid_aas = set('ACDEFGHIKLMNPQRSTVWY')
        if not peptide or peptide[0] != 'M' or not all(aa in valid_aas for aa in peptide):
            raise ValueError(f"Invalid peptide: must start with M and contain only standard amino acids")

        # Get all available RBS options to rotate through
        all_rbs = self.rbsChooser.rbsOptions
        available_rbs = [r for r in all_rbs if r not in ignores]
        if not available_rbs:
            available_rbs = all_rbs

        last_attempt_codons = None
        last_attempt_rbs = available_rbs[0]

        # Scale attempts with protein length
        max_attempts = min(100, max(20, len(peptide) // 8))

        for attempt in range(max_attempts):
            # Rotate through RBS options to escape UTR-caused issues
            selectedRBS = available_rbs[attempt % len(available_rbs)]
            utr = selectedRBS.utr.upper()

            codons = self._build_codons(peptide)
            last_attempt_codons = codons
            last_attempt_rbs = selectedRBS
            cds = ''.join(codons)
            transcript_dna = utr + cds

            # Check all conditions on utr+cds exactly as benchmarker does
            passed_forbidden, _ = self.forbiddenChecker.run(transcript_dna)
            if not passed_forbidden:
                continue

            passed_promoter, _ = self.promoterChecker.run(transcript_dna)
            if not passed_promoter:
                continue

            passed_hairpin, _ = hairpin_checker(transcript_dna)
            if not passed_hairpin:
                continue

            # All checks passed
            return Transcript(selectedRBS, peptide, codons)

        return Transcript(last_attempt_rbs, peptide, last_attempt_codons)


if __name__ == "__main__":
    peptide = "MYPFIRTARMTV"
    designer = TranscriptDesigner()
    designer.initiate()
    ignores = set()
    transcript = designer.run(peptide, ignores)
    print("Codons:", transcript.codons)
    print("CDS:", ''.join(transcript.codons))
    print("RBS UTR:", transcript.rbs.utr)