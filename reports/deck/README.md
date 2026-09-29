# 사내 검토용 슬라이드: 핵심 덱 v1.2와 결과보고 덱 v1.0.1(기술 백업)

- **핵심 덱 v1.2**: `dupilumab_AUCinf_core_deck_v1.2.pptx`와 같은 이름의 PDF. "AUCinf는 1차 평가변수로 적절하지 않다"를 정량 근거로 보이는 짧은 덱(본문 8장(표지 포함) + 별첨 A1~A14 + 보강 B1·B2; 주 모델 = 2020 모델).
  지시 2026-09-29(v1.1 구성, S9 재설계: 1종·2종 오류, v1.2 본문 8장), 사전 등록 `config/prereg_20260929.yaml`·`config/prereg_20260929_oc.yaml`, SPEC §16~17, DECISIONS D-063~D-065. 작성 규칙 `AUTHORING_CORE.md`.
  v1.1 파일(본문 10장)은 v1.2로 대체해 지웠다(커밋 8e5fd7a에 남아 있음).
- **결과보고 덱 v1.0.1(기술 백업)**: `dupilumab_endpoint_results_v1.0.1.pptx`와 PDF(32장). 핵심 덱의 세부 근거를 담은 원본이며 그대로 둔다(파일·추적표·출처 해시 불변).
- 청중: 임상약리·임상개발 담당자. 둘 다 제출 문서가 아니다.

## 핵심 덱 만들기
```
Rscript scripts/63_core_deck_inputs.R          # results/core_deck/(케이스별 창 포착률, 대표 대상자 등; 로컬 .rds 필요, 커밋된 결과와 대조)
Rscript scripts/64_oc_curves.R k2020 2; Rscript scripts/64_oc_curves.R k2016 2   # S7 판정 성능의 새 시험(경계·동일 제품 외 칸)
Rscript scripts/65_oc_curves_summary.R         # results/oc_curves/(저장된 경계·동일 제품 시험 재판정, 1종·2종 오류, 문구 규칙)
Rscript reports/deck/build_deck.R --deck core  # pptx + core_deck_traceability.csv + core_deck_meta*.csv
python3 reports/deck/render_deck.py reports/deck/dupilumab_AUCinf_core_deck_v1.2.pptx <PNG 폴더> \
        --pdf reports/deck/dupilumab_AUCinf_core_deck_v1.2.pdf
Rscript reports/deck/check_deck.R --deck core  # 검사 9종 + 핵심 덱 설계 한도(본문 3줄, 주 그림 60%, 제목 2줄, 글자 하한)
Rscript scripts/60_regulatory_package.R        # 추적 행을 regulatory/traceability.csv(document deck_core_ko)에 합친다
```

## 결과보고 덱(v1.0.1) 만들기
```
Rscript reports/deck/make_template.R          # 16:9 템플릿(처음 한 번; 결과는 template_16x9.pptx로 커밋됨)
Rscript scripts/62_deck_inputs.R               # 분포 그림 입력(results/deck_inputs/; 로컬 대상자 수준 .rds 필요, 커밋된 요약과 대조)
Rscript reports/deck/build_deck.R              # pptx + 추적표(deck_traceability.csv) + 슬라이드 목록(deck_meta*.csv)
python3 reports/deck/render_deck.py reports/deck/dupilumab_endpoint_results_v1.0.1.pptx <PNG 폴더> \
        --pdf reports/deck/dupilumab_endpoint_results_v1.0.1.pdf          # PDF와 슬라이드별 PNG
Rscript reports/deck/check_deck.R              # 자동 검사(실패 시 0이 아닌 종료 코드)
Rscript scripts/60_regulatory_package.R        # 덱 추적 행을 regulatory/traceability.csv(document deck_ko)에 합치고 해시 목록에 덱 파일을 넣는다
```
영문판: `text/en/`에 같은 키 구조의 문구 파일을 두고 `build_deck.R --lang en`.

## 구성
- `lib/deck_lib.R`: 배치(officer), 문구(yaml, 자리표시자), 수치 함수(`d*`, 출처 기록), 넘침 추정, 그림 테마, 글머리표 후처리.
- `lib/deck_facts_common.R`: 설계 상수(config에서 읽음). `slides/sNN_*.R`: 슬라이드별 코드. `text/ko/*.yaml`: 문구(숫자 없음).
- `check_deck.R`: 숫자 추적, 빈 값, 대시, 허용 별첨(결과보고 덱 A5, 핵심 덱 A11) 밖 아토피, 약어 첫 등장, 글자 크기·요점·표 행 수, 출처 해시, M&S 보고서 추적표·key_numbers_en.md와의 대조, 렌더링 배치.
- 핵심 덱: `slides_core/`(공용 도우미 `c00_helpers.R`), `text/ko_core/`, 그림 `figures_core/`. 덱 종류별 설정은 `lib/deck_lib.R`의 `DECK_PROFILES`. 본문 슬라이드의 별첨 안내 줄은 문구 파일의 `xref`.
- `AUTHORING.md`: 작성 규칙(공통). `AUTHORING_CORE.md`: 핵심 덱 추가 규칙.

## 실행 환경(렌더링)
- R 패키지: officer 0.6.4, flextable 0.9.4, ragg 1.2.7, systemfonts 1.0.5와 의존 패키지(renv.lock에 기록).
- 글꼴: Pretendard 1.3.9(SIL Open Font License), https://github.com/orioncactus/pretendard/releases/download/v1.3.9/Pretendard-1.3.9.zip
  (zip SHA-256 04be351a74d6bf7d60c480a3087e51d185485d35a52023142af1df19eb8c428a; `public/static/alternative/*.ttf`를 사용자 글꼴 폴더에 설치).
  pptx는 글꼴 이름만 담는다. Pretendard가 없는 PC의 PowerPoint는 다른 글꼴로 바꿔 보여 준다(PDF에는 글꼴이 들어 있다).
- LibreOffice(Impress, python3-uno)와 PyMuPDF 1.24.10: PDF와 PNG. PDF는 LibreOffice의 한글·영문 사이 자동 간격을 끄고 만든다(PowerPoint 화면과 같게; pptx는 그대로).
