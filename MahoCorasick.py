import ahocorasick
import csv
from pathlib import Path
from typing import Dict, NamedTuple, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class Occurrence(NamedTuple):
    start_index: int
    end_index: int


def get_genome_seqs_from_folder(folder: Path) -> Dict[str, SeqRecord]:
    all_sequences = {}
    for file in folder.glob('*.fasta'):
        with file.open() as opened_file:
            all_sequences[file.stem] = SeqIO.read(opened_file, 'fasta')
    return all_sequences


def get_genome_sequences_matches(genomes: Dict[str, SeqRecord], TEs_to_match: List[SeqRecord]) -> Dict[str, Dict[str, List[Occurrence]]]:
    automaton = ahocorasick.Automaton()
    results: Dict[str, Dict[str, List[Occurrence]]] = {}
    for TE in TEs_to_match:
        automaton.add_word(str(TE.seq), (TE.id, str(TE.seq)))

    for genome_name, genome_seq in genomes.items():

        automaton.make_automaton()

        TE_name_to_occurrences = {}
        for end_index, (TE_name, TE_sequence) in automaton.iter(str(genome_seq.seq)):
            start_index = end_index - len(TE_sequence) + 1

            TE_name_to_occurrences.setdefault(TE_name, []).append(Occurrence(start_index, end_index))

        results[genome_name] = TE_name_to_occurrences

    return results


def export_as_individual_occurrences(matches_and_occurrences: Dict[str, Dict[str, List[Occurrence]]],
                                     output_folder: Path):
    for genome_name, matches in matches_and_occurrences.items():
        with open(output_folder.joinpath(genome_name), mode='w', newline='') as file:
            writer = csv.writer(file)
            for match, occurrences in matches.items():
                for occurrence in occurrences:
                    writer.writerow([match, occurrence.start_index, occurrence.end_index])


def export_as_genome_name_and_match_occurrence_count(matches_and_occurrences: Dict[str, Dict[str, List[Occurrence]]],
                                                     all_TEs: List[str], output_file: Path):
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Fasta Name \ words"] + all_TEs)

        for genome_name, matches in matches_and_occurrences.items():
            for match, occurrences in matches.items():
                writer.writerow([genome_name, match, len(occurrences)])


all_genomes = get_genome_seqs_from_folder(Path(r"C:\Users\Pc\Documents\Maor\ClariTE\test_mahocorasick"))
all_TEs_to_search = SeqIO.parse(r"C:\Users\Pc\Documents\Maor\ClariTE\test_mahocorasick\tes.fna", 'fasta')

all_matches = get_genome_sequences_matches(all_genomes, all_TEs_to_search)

export_as_individual_occurrences(all_matches, Path(r"C:\Users\Pc\Documents\Maor\ClariTE\test_mahocorasick\output"))
