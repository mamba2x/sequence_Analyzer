from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib import parse, request
from urllib.error import HTTPError, URLError

from flask import Flask, render_template, request as flask_request


biology_portal = Flask(__name__)
biology_portal.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024


codon_directory = {
    "UUU": "Phe",
    "UUC": "Phe",
    "UUA": "Leu",
    "UUG": "Leu",
    "UCU": "Ser",
    "UCC": "Ser",
    "UCA": "Ser",
    "UCG": "Ser",
    "UAU": "Tyr",
    "UAC": "Tyr",
    "UAA": "Stop",
    "UAG": "Stop",
    "UGU": "Cys",
    "UGC": "Cys",
    "UGA": "Stop",
    "UGG": "Trp",
    "CUU": "Leu",
    "CUC": "Leu",
    "CUA": "Leu",
    "CUG": "Leu",
    "CCU": "Pro",
    "CCC": "Pro",
    "CCA": "Pro",
    "CCG": "Pro",
    "CAU": "His",
    "CAC": "His",
    "CAA": "Gln",
    "CAG": "Gln",
    "CGU": "Arg",
    "CGC": "Arg",
    "CGA": "Arg",
    "CGG": "Arg",
    "AUU": "Ile",
    "AUC": "Ile",
    "AUA": "Ile",
    "AUG": "Met",
    "ACU": "Thr",
    "ACC": "Thr",
    "ACA": "Thr",
    "ACG": "Thr",
    "AAU": "Asn",
    "AAC": "Asn",
    "AAA": "Lys",
    "AAG": "Lys",
    "AGU": "Ser",
    "AGC": "Ser",
    "AGA": "Arg",
    "AGG": "Arg",
    "GUU": "Val",
    "GUC": "Val",
    "GUA": "Val",
    "GUG": "Val",
    "GCU": "Ala",
    "GCC": "Ala",
    "GCA": "Ala",
    "GCG": "Ala",
    "GAU": "Asp",
    "GAC": "Asp",
    "GAA": "Glu",
    "GAG": "Glu",
    "GGU": "Gly",
    "GGC": "Gly",
    "GGA": "Gly",
    "GGG": "Gly",
}


amino_library = {
    "Ala": {"full_name": "Alanine", "one_letter": "A", "mass": 89.09, "group": "hydrophobic"},
    "Arg": {"full_name": "Arginine", "one_letter": "R", "mass": 174.20, "group": "positively charged"},
    "Asn": {"full_name": "Asparagine", "one_letter": "N", "mass": 132.12, "group": "polar"},
    "Asp": {"full_name": "Aspartic acid", "one_letter": "D", "mass": 133.10, "group": "negatively charged"},
    "Cys": {"full_name": "Cysteine", "one_letter": "C", "mass": 121.16, "group": "polar"},
    "Gln": {"full_name": "Glutamine", "one_letter": "Q", "mass": 146.15, "group": "polar"},
    "Glu": {"full_name": "Glutamic acid", "one_letter": "E", "mass": 147.13, "group": "negatively charged"},
    "Gly": {"full_name": "Glycine", "one_letter": "G", "mass": 75.07, "group": "special"},
    "His": {"full_name": "Histidine", "one_letter": "H", "mass": 155.16, "group": "positively charged"},
    "Ile": {"full_name": "Isoleucine", "one_letter": "I", "mass": 131.18, "group": "hydrophobic"},
    "Leu": {"full_name": "Leucine", "one_letter": "L", "mass": 131.18, "group": "hydrophobic"},
    "Lys": {"full_name": "Lysine", "one_letter": "K", "mass": 146.19, "group": "positively charged"},
    "Met": {"full_name": "Methionine", "one_letter": "M", "mass": 149.21, "group": "hydrophobic"},
    "Phe": {"full_name": "Phenylalanine", "one_letter": "F", "mass": 165.19, "group": "hydrophobic"},
    "Pro": {"full_name": "Proline", "one_letter": "P", "mass": 115.13, "group": "special"},
    "Ser": {"full_name": "Serine", "one_letter": "S", "mass": 105.09, "group": "polar"},
    "Thr": {"full_name": "Threonine", "one_letter": "T", "mass": 119.12, "group": "polar"},
    "Trp": {"full_name": "Tryptophan", "one_letter": "W", "mass": 204.23, "group": "hydrophobic"},
    "Tyr": {"full_name": "Tyrosine", "one_letter": "Y", "mass": 181.19, "group": "polar"},
    "Val": {"full_name": "Valine", "one_letter": "V", "mass": 117.15, "group": "hydrophobic"},
}


@dataclass
class AnalysisSnapshot:
    input_origin: str
    raw_sequence: str
    cleaned_sequence: str
    sequence_family: str
    strand_option: str | None
    detection_story: str
    transcription_story: str
    translation_story: str
    amino_story: str
    protein_story: str
    mrna_result: str
    codon_cards: list[dict[str, Any]]
    residue_cards: list[dict[str, Any]]
    protein_panel: dict[str, Any]
    lookup_cards: list[dict[str, str]]
    warning_text: str | None = None
    error_text: str | None = None


def extract_sequence_payload(uploaded_file: Any, typed_text: str) -> tuple[str, str, str]:
    typed_candidate = typed_text.strip()
    if uploaded_file and uploaded_file.filename:
        file_bytes = uploaded_file.read()
        decoded_payload = file_bytes.decode("utf-8", errors="ignore")
        return decoded_payload, uploaded_file.filename, "uploaded file"
    if typed_candidate:
        return typed_candidate, "typed_input", "text area"
    return "", "", ""


def normalize_sequence_text(source_blob: str) -> str:
    normalized_lines = []
    for source_line in source_blob.splitlines():
        trimmed_line = source_line.strip()
        if not trimmed_line or trimmed_line.startswith(">"):
            continue
        normalized_lines.append(trimmed_line)
    merged_fragment = "".join(normalized_lines)
    return "".join(character_unit for character_unit in merged_fragment.upper() if not character_unit.isspace())


def identify_sequence_family(sequence_probe: str) -> tuple[str, str]:
    if not sequence_probe:
        return "invalid", "The application could not find any letters that look like a sequence after cleaning the input."

    dna_symbols = set("ACGT")
    rna_symbols = set("ACGU")
    symbol_inventory = set(sequence_probe)
    dna_match = symbol_inventory.issubset(dna_symbols)
    rna_match = symbol_inventory.issubset(rna_symbols)

    if dna_match and "U" not in symbol_inventory:
        return "DNA", (
            "This was identified as DNA because every letter belongs to the DNA alphabet "
            "(A, C, G, and T), and the sequence uses T instead of U."
        )
    if rna_match and "U" in symbol_inventory and "T" not in symbol_inventory:
        return "RNA", (
            "This was identified as RNA because every letter belongs to the RNA alphabet "
            "(A, C, G, and U), and the sequence uses U instead of T."
        )

    invalid_symbols = sorted(symbol_inventory.difference(set("ACGTU")))
    if invalid_symbols:
        joined_faults = ", ".join(invalid_symbols)
        return "invalid", (
            f"This input is invalid because it contains letters that do not belong to DNA or RNA: {joined_faults}."
        )

    return "invalid", (
        "This input mixes T and U together, so it does not cleanly match DNA or RNA."
    )


def transcribe_sequence(sequence_family: str, sequence_probe: str, strand_option: str | None) -> tuple[str, str]:
    if sequence_family == "RNA":
        transcript_result = sequence_probe
        transcript_story = (
            "Transcription is the step where cells make mRNA from DNA. In this case the input was already RNA, "
            "so the application treated it as the messenger RNA sequence directly instead of converting it again."
        )
        return transcript_result, transcript_story

    if strand_option == "template":
        transcription_map = str.maketrans({"A": "U", "T": "A", "C": "G", "G": "C"})
        transcript_result = sequence_probe.translate(transcription_map)
        transcript_story = (
            "Transcription turns DNA information into mRNA. Because you marked the input as the template strand, "
            "the program built a complementary RNA copy from your exact letters: A became U, T became A, C became G, and G became C."
        )
        return transcript_result, transcript_story

    transcript_result = sequence_probe.replace("T", "U")
    transcript_story = (
        "Transcription turns DNA information into mRNA. Because you marked the input as the coding "
        "strand, the program kept the same order of letters and changed each T to U to match RNA."
    )
    return transcript_result, transcript_story


def translate_mrna(transcript_result: str) -> tuple[list[dict[str, Any]], list[str], str]:
    opening_index = transcript_result.find("AUG")
    codon_cards: list[dict[str, Any]] = []
    residue_chain: list[str] = []

    if opening_index == -1:
        translation_story = (
            "Translation is the step where the cell reads mRNA three letters at a time to build amino acids. "
            "No start codon (AUG) was found here, so the program could not begin a protein-coding translation."
        )
        return codon_cards, residue_chain, translation_story

    trailing_note = ""
    for codon_cursor in range(opening_index, len(transcript_result), 3):
        codon_triplet = transcript_result[codon_cursor : codon_cursor + 3]
        if len(codon_triplet) < 3:
            trailing_note = (
                f" The final leftover bases '{codon_triplet}' were ignored because translation only reads full groups of three."
            )
            break

        amino_code = codon_directory.get(codon_triplet, "Unknown")
        if amino_code == "Stop":
            codon_cards.append(
                {
                    "position": len(codon_cards) + 1,
                    "codon": codon_triplet,
                    "amino_short": "Stop",
                    "amino_name": "Termination signal",
                    "note": "This codon tells translation to stop.",
                }
            )
            translation_story = (
                "Translation reads the mRNA in codons, which are groups of three bases. The program started at the first AUG "
                "start codon, converted each codon into its amino acid, and stopped when it reached a stop codon."
                + trailing_note
            )
            return codon_cards, residue_chain, translation_story

        amino_entry = amino_library[amino_code]
        residue_chain.append(amino_code)
        codon_cards.append(
            {
                "position": len(codon_cards) + 1,
                "codon": codon_triplet,
                "amino_short": amino_code,
                "amino_name": amino_entry["full_name"],
                "note": "This codon adds one amino acid to the growing chain.",
            }
        )

    translation_story = (
        "Translation reads the mRNA in codons, which are groups of three bases. The program started at the first AUG start codon "
        "and converted each full codon into an amino acid. No stop codon appeared before the sequence ended, so the chain is shown as incomplete."
        + trailing_note
    )
    return codon_cards, residue_chain, translation_story


def build_residue_cards(residue_chain: list[str]) -> tuple[list[dict[str, str]], str]:
    residue_cards = []
    for residue_index, residue_code in enumerate(residue_chain, start=1):
        residue_record = amino_library[residue_code]
        residue_cards.append(
            {
                "position": str(residue_index),
                "full_name": residue_record["full_name"],
                "three_letter": residue_code,
                "one_letter": residue_record["one_letter"],
                "group": residue_record["group"],
            }
        )

    if residue_cards:
        amino_story = (
            "Amino acids are the small chemical building blocks that get linked together during translation. "
            "This ordered list is the polypeptide chain, which is the direct starting material from which a protein can form."
        )
    else:
        amino_story = (
            "Amino acids are the building blocks of proteins, but none were produced here because translation never started or did not yield a codon-to-amino-acid chain."
        )
    return residue_cards, amino_story


def characterize_protein(residue_chain: list[str]) -> dict[str, Any]:
    if not residue_chain:
        return {
            "protein_sequence": "",
            "residue_total": 0,
            "mass_estimate": 0.0,
            "composition_text": "No protein-like chain could be characterized because no amino acids were produced.",
            "property_table": [],
        }

    one_letter_sequence = "".join(amino_library[residue_code]["one_letter"] for residue_code in residue_chain)
    property_counter = Counter(amino_library[residue_code]["group"] for residue_code in residue_chain)
    property_table = [
        {"label": group_name.title(), "count": property_counter.get(group_name, 0)}
        for group_name in [
            "hydrophobic",
            "polar",
            "positively charged",
            "negatively charged",
            "special",
        ]
    ]
    composition_parts = [f"{row_entry['count']} {row_entry['label'].lower()} residues" for row_entry in property_table]
    return {
        "protein_sequence": one_letter_sequence,
        "residue_total": len(residue_chain),
        "mass_estimate": round(sum(amino_library[residue_code]["mass"] for residue_code in residue_chain), 2),
        "composition_text": ", ".join(composition_parts),
        "property_table": property_table,
    }


def search_uniprot(one_letter_sequence: str) -> tuple[list[dict[str, str]], str]:
    if not one_letter_sequence:
        return [], (
            "Proteins are folded chains of amino acids. Because no amino-acid chain was produced, there was no protein sequence to send to UniProt for matching."
        )

    query_phrase = (
        f"https://rest.uniprot.org/uniprotkb/search?query=sequence:{parse.quote(one_letter_sequence)}"
        "&fields=protein_name,organism_name,cc_function&format=json&size=5"
    )

    try:
        with request.urlopen(query_phrase, timeout=10) as lookup_response:
            payload_text = lookup_response.read().decode("utf-8")
    except (HTTPError, URLError, TimeoutError):
        return [], (
            "Proteins are functional molecules made from amino-acid chains. The application prepared a UniProt lookup "
            "for your sequence, but the database could not be reached from this environment, so no live matches are shown."
        )

    import json

    payload_block = json.loads(payload_text)
    result_cards = []
    for hit_entry in payload_block.get("results", []):
        protein_block = hit_entry.get("proteinDescription", {})
        recommended_block = protein_block.get("recommendedName", {})
        full_name_block = recommended_block.get("fullName", {})
        comment_stack = hit_entry.get("comments", [])
        function_blurb = "Function information was not returned."
        for comment_unit in comment_stack:
            if comment_unit.get("commentType") == "FUNCTION":
                comment_texts = comment_unit.get("texts", [])
                if comment_texts:
                    function_blurb = comment_texts[0].get("value", function_blurb)
                    break

        result_cards.append(
            {
                "protein_name": full_name_block.get("value", "Unnamed protein"),
                "organism_name": hit_entry.get("organism", {}).get("scientificName", "Unknown organism"),
                "function_text": function_blurb,
                "accession_id": hit_entry.get("primaryAccession", "N/A"),
            }
        )

    if result_cards:
        protein_story = (
            "A protein is the finished molecule that can form when a polypeptide chain folds into a working shape. "
            "The app summarized your translated chain and then searched UniProt for real biological records whose sequence matches this protein."
        )
    else:
        protein_story = (
            "A protein is the working molecule built from a polypeptide chain. The app searched UniProt with the translated sequence, "
            "but no direct matches were returned for the current query."
        )
    return result_cards, protein_story


def analyze_submission(raw_blob: str, input_origin: str, strand_option: str | None) -> AnalysisSnapshot:
    cleaned_sequence = normalize_sequence_text(raw_blob)
    sequence_family, detection_story = identify_sequence_family(cleaned_sequence)

    if sequence_family == "invalid":
        return AnalysisSnapshot(
            input_origin=input_origin,
            raw_sequence=raw_blob,
            cleaned_sequence=cleaned_sequence,
            sequence_family=sequence_family,
            strand_option=strand_option,
            detection_story=detection_story,
            transcription_story="Transcription was not run because the input did not qualify as valid DNA or RNA.",
            translation_story="Translation was not run because there was no valid mRNA sequence to read.",
            amino_story="Amino-acid reporting was skipped because translation did not produce a chain.",
            protein_story="Protein characterization was skipped because there was no translated chain to analyze.",
            mrna_result="",
            codon_cards=[],
            residue_cards=[],
            protein_panel=characterize_protein([]),
            lookup_cards=[],
            error_text="Please provide a valid DNA or RNA sequence using A, C, G, T, or U.",
        )

    transcript_result, transcription_story = transcribe_sequence(sequence_family, cleaned_sequence, strand_option)
    codon_cards, residue_chain, translation_story = translate_mrna(transcript_result)
    residue_cards, amino_story = build_residue_cards(residue_chain)
    protein_panel = characterize_protein(residue_chain)
    lookup_cards, protein_story = search_uniprot(protein_panel["protein_sequence"])

    warning_text = None
    if sequence_family == "DNA" and not strand_option:
        warning_text = "No DNA strand type was selected, so the program assumed the coding (non-template) strand."

    return AnalysisSnapshot(
        input_origin=input_origin,
        raw_sequence=raw_blob,
        cleaned_sequence=cleaned_sequence,
        sequence_family=sequence_family,
        strand_option=strand_option or "coding",
        detection_story=detection_story,
        transcription_story=transcription_story,
        translation_story=translation_story,
        amino_story=amino_story,
        protein_story=protein_story,
        mrna_result=transcript_result,
        codon_cards=codon_cards,
        residue_cards=residue_cards,
        protein_panel=protein_panel,
        lookup_cards=lookup_cards,
        warning_text=warning_text,
    )


@biology_portal.route("/", methods=["GET", "POST"])
def home_page() -> str:
    rendered_snapshot = None

    if flask_request.method == "POST":
        uploaded_file = flask_request.files.get("sequence_file")
        typed_text = flask_request.form.get("sequence_text", "")
        strand_option = flask_request.form.get("dna_strand")
        raw_blob, input_label, input_origin = extract_sequence_payload(uploaded_file, typed_text)

        if raw_blob:
            rendered_snapshot = analyze_submission(raw_blob, input_origin or input_label, strand_option)
        else:
            rendered_snapshot = AnalysisSnapshot(
                input_origin="none",
                raw_sequence="",
                cleaned_sequence="",
                sequence_family="invalid",
                strand_option=None,
                detection_story="No sequence was submitted yet.",
                transcription_story="",
                translation_story="",
                amino_story="",
                protein_story="",
                mrna_result="",
                codon_cards=[],
                residue_cards=[],
                protein_panel=characterize_protein([]),
                lookup_cards=[],
                error_text="Enter a sequence in the text box or upload a text/FASTA file first.",
            )

    return render_template("index.html", rendered_snapshot=rendered_snapshot)


if __name__ == "__main__":
    template_folder = Path(__file__).with_name("templates")
    static_folder = Path(__file__).with_name("static")
    if not template_folder.exists() or not static_folder.exists():
        raise SystemExit("Required folders are missing. Expected 'templates' and 'static' next to app.py.")
    biology_portal.run(host="127.0.0.1", port=5001, debug=True)

