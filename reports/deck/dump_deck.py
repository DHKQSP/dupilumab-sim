#!/usr/bin/env python3
"""dump_deck.py: pptx의 슬라이드별 글자(도형 이름별)와 발표자 노트를 텍스트로 뽑는다(검토용).

사용법: python3 reports/deck/dump_deck.py deck.pptx OUT.md
"""
import re, sys, zipfile
from lxml import etree

NS = {"a": "http://schemas.openxmlformats.org/drawingml/2006/main", "p": "http://schemas.openxmlformats.org/presentationml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def para_text(node):
    out = []
    for p in node.iterfind(".//a:p", NS):
        t = "".join(x.text or "" for x in p.iterfind(".//a:t", NS))
        if t.strip():
            out.append(t)
    return out


def main(pptx, out):
    z = zipfile.ZipFile(pptx)
    pres = etree.fromstring(z.read("ppt/presentation.xml")); rels = etree.fromstring(z.read("ppt/_rels/presentation.xml.rels"))
    rmap = {r.get("Id"): r.get("Target") for r in rels}
    ids = [s.get("{%s}id" % NS["r"]) for s in pres.iterfind(".//p:sldIdLst/p:sldId", NS)]
    lines = []
    for n, rid in enumerate(ids, start=1):
        path = "ppt/" + rmap[rid]; x = etree.fromstring(z.read(path))
        lines.append(f"\n## Slide {n} ({path})")
        for sh in x.iterfind(".//p:spTree/*", NS):
            name_el = sh.find(".//p:cNvPr", NS)
            name = name_el.get("name") if name_el is not None else "?"
            txt = para_text(sh)
            if name.startswith("figure_"):
                lines.append(f"- [{name}] (image)")
            elif txt:
                if sh.tag.endswith("graphicFrame"):
                    rows = []
                    for tr in sh.iterfind(".//a:tr", NS):
                        rows.append(" | ".join(" ".join(para_text(tc)) for tc in tr.iterfind("a:tc", NS)))
                    lines.append(f"- [{name}] table:\n    " + "\n    ".join(rows))
                else:
                    lines.append(f"- [{name}] " + " / ".join(txt))
        rp = re.sub(r"slides/", "slides/_rels/", path) + ".rels"
        if rp in z.namelist():
            rr = etree.fromstring(z.read(rp))
            for r in rr:
                if "notesSlide" in r.get("Target"):
                    npath = "ppt/" + r.get("Target").replace("../", "")
                    nx = etree.fromstring(z.read(npath))
                    body = [t for sp in nx.iterfind(".//p:sp", NS) if sp.find(".//p:ph[@type='body']", NS) is not None for t in para_text(sp)]
                    if body:
                        lines.append("- [notes] " + " / ".join(body))
    open(out, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    print(f"{len(ids)} slides -> {out}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
