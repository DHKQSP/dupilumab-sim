#!/usr/bin/env python3
"""check_render.py: 렌더링된 PDF에서 글자 넘침·겹침을 찾는다(넘침 추정만으로는 놓치는 것을 실제 그림으로 확인).

pptx의 도형 위치(글상자, 표, 그림)와 PDF의 글줄 위치(PyMuPDF)를 맞대어 본다.
- overflow: 글줄이 자기 도형 상자 밖으로 2 pt 넘게 나감(자기 도형 = 그 글줄의 글자를 담은 도형)
- overlap_text: 서로 다른 도형의 글줄끼리 겹침
- overlap_figure: 글줄이 그림 상자와 겹침(그림 위에 얹은 글상자)
- margin: 글줄이 슬라이드 가장자리 0.3 in 안쪽까지 들어옴
- unmatched: 어느 도형에도 대응하지 않는 글줄(검사 불가, 보고만 함)
그림 안 글자(ragg PNG)는 그림 자체라 여기서 보지 않는다(육안 확인).

사용법: python3 reports/deck/check_render.py deck.pptx deck.pdf OUT.csv
종료 코드: overflow·overlap·margin이 하나라도 있으면 1.
"""
import csv, re, sys, zipfile
from lxml import etree

NS = {"a": "http://schemas.openxmlformats.org/drawingml/2006/main", "p": "http://schemas.openxmlformats.org/presentationml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
EMU_PT = 12700.0
TOL = 2.0          # pt
MARGIN = 0.3 * 72  # pt


def norm(s):
    s = s.replace(chr(0x2011), "-")
    s = re.sub("[" + chr(0x200b) + chr(0x2060) + chr(0x2063) + chr(0xfeff) + "]", "", s)
    s = re.sub("^[\\s" + chr(0x25cf) + chr(0x2022) + chr(0x2013) + "]+", "", s)     # 글머리표 문자
    return re.sub(r"\s+", "", s)


def para_texts(node):
    out = []
    for p in node.iterfind(".//a:p", NS):
        t = "".join(x.text or "" for x in p.iterfind(".//a:t", NS))
        t = re.sub(chr(0x2063) + "B[0-9]" + chr(0x2063), "", t)
        if t.strip():
            out.append(t)
    return out


def shapes_of(z, path):
    x = etree.fromstring(z.read(path)); out = []
    for sh in x.iterfind(".//p:spTree/*", NS):
        tag = etree.QName(sh).localname
        if tag not in ("sp", "graphicFrame", "pic"):
            continue
        name_el = sh.find(".//p:cNvPr", NS); name = name_el.get("name") if name_el is not None else "?"
        xf = sh.find("p:spPr/a:xfrm", NS) if tag in ("sp", "pic") else sh.find("p:xfrm", NS)
        if xf is None:
            continue
        off = xf.find("a:off", NS); ext = xf.find("a:ext", NS)
        box = (int(off.get("x")) / EMU_PT, int(off.get("y")) / EMU_PT, (int(off.get("x")) + int(ext.get("cx"))) / EMU_PT, (int(off.get("y")) + int(ext.get("cy"))) / EMU_PT)
        kind = "figure" if tag == "pic" else ("table" if tag == "graphicFrame" else "text")
        paras = para_texts(sh) if kind != "figure" else []
        if kind == "table":   # 셀 단위 문단
            paras = [" ".join(para_texts(tc)) for tc in sh.iterfind(".//a:tc", NS)]
        full = norm("".join(paras))
        out.append(dict(name=name, kind=kind, box=box, text=full, paras=[norm(p) for p in paras if norm(p)]))
    return out


def inter(a, b):
    w = min(a[2], b[2]) - max(a[0], b[0]); h = min(a[3], b[3]) - max(a[1], b[1])
    return (w, h) if w > 0 and h > 0 else (0, 0)


def main(pptx, pdf, out):
    import fitz
    z = zipfile.ZipFile(pptx)
    pres = etree.fromstring(z.read("ppt/presentation.xml")); rels = etree.fromstring(z.read("ppt/_rels/presentation.xml.rels"))
    rmap = {r.get("Id"): r.get("Target") for r in rels}
    ids = [s.get("{%s}id" % NS["r"]) for s in pres.iterfind(".//p:sldIdLst/p:sldId", NS)]
    sz = pres.find("p:sldSz", NS); W = int(sz.get("cx")) / EMU_PT; H = int(sz.get("cy")) / EMU_PT
    doc = fitz.open(pdf)
    if doc.page_count != len(ids):
        raise SystemExit(f"page count {doc.page_count} != slides {len(ids)}")
    rows = []
    for n, rid in enumerate(ids, start=1):
        shapes = shapes_of(z, "ppt/" + rmap[rid]); page = doc[n - 1]
        sx = page.rect.width / W; sy = page.rect.height / H
        lines = []
        for b in page.get_text("dict")["blocks"]:
            for ln in b.get("lines", []):
                t = "".join(s["text"] for s in ln["spans"])
                if not norm(t):
                    continue
                x0, y0, x1, y1 = ln["bbox"]; lines.append(dict(text=t, nt=norm(t), box=(x0 / sx, y0 / sy, x1 / sx, y1 / sy)))
        # 글줄마다 자기 도형 찾기: 글자가 들어 있는 도형 중 상자와 가장 많이 겹치는 것(없으면 가장 가까운 것)
        for ln in lines:
            cands = [s for s in shapes if s["kind"] != "figure" and ln["nt"] in s["text"]]
            if not cands:
                rows.append(dict(slide=n, kind="unmatched", shape="", text=ln["text"][:60], detail="line not found in any shape")); ln["owner"] = None; continue

            def score(s):
                w, h = inter(ln["box"], s["box"]); area = w * h
                cx = (ln["box"][0] + ln["box"][2]) / 2; cy = (ln["box"][1] + ln["box"][3]) / 2
                d = max(0, s["box"][0] - cx, cx - s["box"][2]) + max(0, s["box"][1] - cy, cy - s["box"][3])
                return (area, -d)
            own = max(cands, key=score); ln["owner"] = own["name"]
            b, s = ln["box"], own["box"]
            ex = max(s[0] - b[0], b[2] - s[2], s[1] - b[1], b[3] - s[3])
            if ex > TOL:
                side = [k for k, v in (("left", s[0] - b[0]), ("right", b[2] - s[2]), ("top", s[1] - b[1]), ("bottom", b[3] - s[3])) if v > TOL]
                rows.append(dict(slide=n, kind="overflow", shape=own["name"], text=ln["text"][:60], detail=f"{'/'.join(side)} by {ex:.1f} pt"))
            mg = 0.1 * 72 if own["name"].startswith("footer") else MARGIN   # 바닥글은 가장자리 0.1 in까지 허용
            if b[0] < mg - TOL or b[1] < mg - TOL or b[2] > W - mg + TOL or b[3] > H - mg + TOL:
                rows.append(dict(slide=n, kind="margin", shape=own["name"], text=ln["text"][:60], detail="within 0.3 in of the slide edge"))
            for f in shapes:
                if f["kind"] == "figure":
                    w, h = inter(b, f["box"])
                    if w > TOL and h > TOL:
                        rows.append(dict(slide=n, kind="overlap_figure", shape=own["name"], text=ln["text"][:60], detail=f"overlaps {f['name']} ({w:.0f}x{h:.0f} pt)"))
        for i in range(len(lines)):
            for j in range(i + 1, len(lines)):
                a, c = lines[i], lines[j]
                if a.get("owner") and c.get("owner") and a["owner"] != c["owner"]:
                    w, h = inter(a["box"], c["box"])
                    if w > TOL and h > TOL:
                        rows.append(dict(slide=n, kind="overlap_text", shape=f"{a['owner']} x {c['owner']}", text=f"{a['text'][:30]} | {c['text'][:30]}", detail=f"{w:.0f}x{h:.0f} pt"))
    with open(out, "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=["slide", "kind", "shape", "text", "detail"]); wr.writeheader(); wr.writerows(rows)
    bad = [r for r in rows if r["kind"] != "unmatched"]
    print(f"{len(ids)} slides; {len(bad)} layout problems, {len(rows) - len(bad)} unmatched lines -> {out}")
    for r in bad:
        print(f"  slide {r['slide']}: {r['kind']} [{r['shape']}] {r['detail']} :: {r['text']}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3]))
