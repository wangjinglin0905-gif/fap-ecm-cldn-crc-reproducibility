import json
import argparse
import urllib.parse
import urllib.request
from pathlib import Path


PMIDS = [
    "25706627",  # Isella
    "25706628",  # Calon
    "32961620",  # Gonzalez-Gualda
    "38831121",  # SenNet
    "33531494",  # Zimmerman
    "34584091",  # Squair
    "33159064",  # Avila Cobos
    "33168968",  # Armingol
    "22496421",  # Freund
    "41535933",  # Ge 2026
    "42418705",  # Yang 2026
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Verify the frozen bibliography metadata through NCBI E-utilities."
    )
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()
    query = urllib.parse.urlencode(
        {"db": "pubmed", "id": ",".join(PMIDS), "retmode": "json"}
    )
    request = urllib.request.Request(
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{query}",
        headers={"User-Agent": "Codex manuscript citation verification/1.0"},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        payload = json.load(response)

    records = []
    for pmid in payload["result"]["uids"]:
        item = payload["result"][pmid]
        doi = next(
            (x["value"] for x in item.get("articleids", []) if x.get("idtype") == "doi"),
            "",
        )
        records.append(
            {
                "pmid": pmid,
                "title": item.get("title", ""),
                "journal": item.get("fulljournalname", ""),
                "pubdate": item.get("pubdate", ""),
                "volume": item.get("volume", ""),
                "issue": item.get("issue", ""),
                "pages": item.get("pages", ""),
                "doi": doi,
                "authors": [x.get("name", "") for x in item.get("authors", [])],
            }
        )

    target = args.output_json
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8")
    for record in records:
        print(
            f"{record['pmid']} | {record['title']} | {record['journal']} | "
            f"{record['pubdate']} | {record['volume']}({record['issue']}):{record['pages']} | "
            f"{record['doi']} | {'; '.join(record['authors'][:4])}"
        )


if __name__ == "__main__":
    main()
