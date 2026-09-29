#!/usr/bin/env python3
"""render_deck.py: 슬라이드 pptx를 PDF와 슬라이드별 PNG로 변환한다(열람용 PDF, 육안 확인·자동 검사용 PNG).

기본 경로: LibreOffice(headless)를 UNO로 띄워 pptx를 열고, 모든 문단의 '아시아·라틴 문자 사이 간격 자동 조정'
(ParaIsCharacterDistance)을 끈 뒤 PDF로 내보낸다. LibreOffice는 이 설정을 켜 두어 한글과 영문·숫자 사이에 간격을 더
넣어 그리는데, PowerPoint는 그렇게 그리지 않으므로 PDF가 PowerPoint 화면과 같아지게 하려는 것이다(pptx 파일은 바꾸지 않음).
UNO(python3-uno)를 쓸 수 없으면 soffice --convert-to pdf로 바꾼다. PNG는 PyMuPDF로 그린다.
사용자 프로필은 임시 폴더를 쓴다(동시 실행 가능). 환경변수 DECK_SOFFICE_HELPER가 soffice.py 도우미(소켓 제한 샌드박스용
LD_PRELOAD 심)를 가리키면 그 환경을 쓴다.

사용법: python3 reports/deck/render_deck.py deck.pptx OUTDIR [--dpi 110] [--pages 3,5] [--pdf 경로] [--plain]
"""
import argparse, importlib.util, os, subprocess, sys, tempfile, time
from pathlib import Path


def soffice_env():
    helper = os.environ.get("DECK_SOFFICE_HELPER")
    if helper:
        spec = importlib.util.spec_from_file_location("deck_soffice", helper)
        mod = importlib.util.module_from_spec(spec); sys.path.insert(0, str(Path(helper).parent.parent)); spec.loader.exec_module(mod)
        return mod.get_soffice_env(), mod
    return dict(os.environ, SAL_USE_VCLPLUGIN="svp"), None


def pdf_uno(pptx: Path, pdf: Path) -> int:
    import uno
    from com.sun.star.beans import PropertyValue
    env, _ = soffice_env()
    prof = tempfile.mkdtemp(prefix="lo_uno_"); port = 20000 + (os.getpid() * 7) % 20000
    p = subprocess.Popen(["soffice", "--headless", "--invisible", "--norestore", "--nologo", f"-env:UserInstallation=file://{prof}",
                          f"--accept=socket,host=127.0.0.1,port={port};urp;StarOffice.ComponentContext"], env=env,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        local = uno.getComponentContext()
        resolver = local.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local)
        ctx = None
        for _ in range(120):
            try:
                ctx = resolver.resolve(f"uno:socket,host=127.0.0.1,port={port};urp;StarOffice.ComponentContext"); break
            except Exception:
                time.sleep(0.5)
        if ctx is None:
            raise RuntimeError("could not connect to LibreOffice")
        desktop = ctx.ServiceManager.createInstanceWithContext("com.sun.star.frame.Desktop", ctx)

        def pv(n, v):
            x = PropertyValue(); x.Name = n; x.Value = v; return x
        doc = desktop.loadComponentFromURL(uno.systemPathToFileUrl(str(pptx)), "_blank", 0, (pv("Hidden", True),))
        n = 0

        def fix(text):
            nonlocal n
            en = text.createEnumeration()
            while en.hasMoreElements():
                para = en.nextElement()
                try:
                    para.ParaIsCharacterDistance = False; n += 1
                except Exception:
                    pass
        for i in range(doc.DrawPages.Count):
            page = doc.DrawPages.getByIndex(i)
            for j in range(page.Count):
                sh = page.getByIndex(j)
                if sh.ShapeType == "com.sun.star.drawing.TableShape":
                    m = sh.Model
                    for r in range(m.Rows.Count):
                        for c in range(m.Columns.Count):
                            cell = m.getCellByPosition(c, r)
                            try:
                                fix(cell.Text if hasattr(cell, "Text") else cell)
                            except Exception:
                                fix(cell)
                else:
                    try:
                        fix(sh.Text)
                    except Exception:
                        pass
        doc.storeToURL(uno.systemPathToFileUrl(str(pdf)), (pv("FilterName", "impress_pdf_Export"),))
        doc.close(True)
        try:
            desktop.terminate()
        except Exception:
            pass
        return n
    finally:
        try:
            p.wait(timeout=30)
        except Exception:
            p.kill()


def pdf_plain(pptx: Path, outdir: Path) -> Path:
    env, mod = soffice_env()
    with tempfile.TemporaryDirectory(prefix="lo_profile_") as prof:
        argv = ["--headless", f"-env:UserInstallation=file://{prof}", "--convert-to", "pdf", "--outdir", str(outdir), str(pptx)]
        res = mod.run_soffice(argv, capture_output=True, text=True, timeout=600) if mod else subprocess.run(["soffice"] + argv, capture_output=True, text=True, timeout=600, env=env)
    pdf = outdir / (pptx.stem + ".pdf")
    if not pdf.exists():
        raise SystemExit(f"PDF not produced:\n{res.stdout}\n{res.stderr}")
    return pdf


def pdf_png(pdf: Path, outdir: Path, dpi: int, pages=None):
    import fitz  # PyMuPDF
    doc = fitz.open(pdf); width = len(str(doc.page_count)); out = []
    for i, page in enumerate(doc, start=1):
        if pages and i not in pages:
            continue
        f = outdir / f"slide-{i:0{width}d}.png"; page.get_pixmap(dpi=dpi).save(f); out.append(f)
    return doc.page_count, out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("pptx"); ap.add_argument("outdir"); ap.add_argument("--dpi", type=int, default=110); ap.add_argument("--pages", default="")
    ap.add_argument("--pdf", default=""); ap.add_argument("--plain", action="store_true")
    a = ap.parse_args()
    pptx = Path(a.pptx).resolve(); outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    pdf = Path(a.pdf).resolve() if a.pdf else (outdir / (pptx.stem + ".pdf"))
    how = "plain"
    if not a.plain:
        try:
            n = pdf_uno(pptx, pdf); how = f"uno ({n} paragraphs without Asian/Latin autospacing)"
        except Exception as e:
            print(f"UNO conversion unavailable ({e}); falling back to plain conversion", file=sys.stderr)
    if how == "plain":
        tmp = pdf_plain(pptx, outdir)
        if tmp != pdf:
            tmp.replace(pdf)
    pages = {int(x) for x in a.pages.split(",") if x.strip()} or None
    n, files = pdf_png(pdf, outdir, a.dpi, pages)
    print(f"pdf: {pdf} ({n} pages; {how})")
    for f in files:
        print(f)
