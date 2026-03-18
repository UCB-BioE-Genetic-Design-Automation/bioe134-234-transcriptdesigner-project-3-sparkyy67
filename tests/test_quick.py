from genedesign.transcript_designer import TranscriptDesigner
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.codon_checker import CodonChecker

def parse_fasta_first_n(fasta_file, n=10):
    sequences = {}
    current_gene = None
    current_sequence = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if len(sequences) >= n:
                break
            line = line.strip()
            if line.startswith(">"):
                if current_gene:
                    sequences[current_gene] = ''.join(current_sequence)
                parts = line.split()
                gene_name = None
                for part in parts:
                    if part.startswith("GN="):
                        gene_name = part.split("=")[1]
                        break
                if not gene_name:
                    gene_name = line.split('|')[2].split(' ')[0]
                current_gene = gene_name
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_gene and len(sequences) < n:
            sequences[current_gene] = ''.join(current_sequence)
    return sequences

if __name__ == "__main__":
    fasta_file = "tests/benchmarking/uniprotkb_proteome_UP000054015_2024_09_24.fasta"
    designer = TranscriptDesigner()
    designer.initiate()
    forbidden_checker = ForbiddenSequenceChecker()
    forbidden_checker.initiate()
    promoter_checker = PromoterChecker()
    promoter_checker.initiate()
    codon_checker = CodonChecker()
    codon_checker.initiate()
    proteome = parse_fasta_first_n(fasta_file, n=10)
    passed = 0
    failed = 0
    for gene, protein in proteome.items():
        try:
            ignores = set()
            transcript = designer.run(protein, ignores)
            cds = ''.join(transcript.codons)
            issues = []
            ok, site = forbidden_checker.run(cds)
            if not ok:
                issues.append(f"Forbidden: {site}")
            ok, _ = hairpin_checker(cds)
            if not ok:
                issues.append("Hairpin")
            ok, promoter = promoter_checker.run(cds)
            if not ok:
                issues.append(f"Promoter")
            ok, div, rare, cai = codon_checker.run(transcript.codons)
            if not ok:
                issues.append(f"Codon: div={div:.2f} rare={rare} cai={cai:.2f}")
            if issues:
                print(f"FAIL {gene}: {', '.join(issues)}")
                failed += 1
            else:
                print(f"PASS {gene}")
                passed += 1
        except Exception as e:
            print(f"ERROR {gene}: {e}")
            failed += 1
    print(f"\n{passed} passed, {failed} failed out of {len(proteome)} proteins")