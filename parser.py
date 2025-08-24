import re
from pypdf import PdfReader


def get_cas_numbers_without_braces(text):
    "matches all CAS of form [12-34-56] and returns list with CAS of form ['12-34-56', ... ]"

    res = []
    matches = re.findall(cas_pattern, text)
    for cas in matches:
        res.append(cas.replace("[", "").replace("]", ""))
    return res


if __name__ == "__main__":
    cas_pattern = r"\d{2,7}-\d{2}-\d"

    reader = PdfReader("input.pdf")
    seen = set()

    # only those pages which contain general chemicals for a particular pdf
    pages_with_general_chemicals = reader.pages[18:76]

    for page in pages_with_general_chemicals:
        found_cas_numbers = get_cas_numbers_without_braces(page.extract_text())

        # add cas to a set to avoid duplicates
        for cas in found_cas_numbers:
            seen.add(cas)

    with open("actual_cas.txt","w") as output:
        for cas in seen:
            output.write(cas + "\n")