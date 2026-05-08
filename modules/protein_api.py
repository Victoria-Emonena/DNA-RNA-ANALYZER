"""
modules/protein_api.py — UniProt REST API integration
"""
import urllib.request
import urllib.parse
import json


UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


def search_uniprot(query: str, max_results: int = 5) -> list[dict]:
    """
    Search UniProt for proteins matching the query string.
    Returns a list of result dicts with name, organism, function, accession.
    """
    params = urllib.parse.urlencode({
        "query": query,
        "format": "json",
        "size": max_results,
        "fields": "accession,protein_name,organism_name,cc_function",
    })
    url = f"{UNIPROT_SEARCH_URL}?{params}"

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
    except Exception as e:
        return [{"error": str(e)}]

    results = []
    for entry in data.get("results", []):
        accession = entry.get("primaryAccession", "N/A")
        # Protein name
        pn = entry.get("proteinDescription", {})
        rec = pn.get("recommendedName", {})
        full_name = rec.get("fullName", {}).get("value", "")
        if not full_name:
            sub = pn.get("submissionNames", [{}])
            full_name = sub[0].get("fullName", {}).get("value", "Unknown") if sub else "Unknown"

        organism = entry.get("organism", {}).get("scientificName", "N/A")

        # Function from comments
        function_text = "N/A"
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    function_text = texts[0].get("value", "N/A")[:300]
                    break

        results.append({
            "accession": accession,
            "name": full_name,
            "organism": organism,
            "function": function_text,
            "url": f"https://www.uniprot.org/uniprotkb/{accession}",
        })

    return results
