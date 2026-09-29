# 결과보고 슬라이드 작성 규칙 (reports/deck)

슬라이드는 R officer로 만든다. 슬라이드 하나 = `slides/sNN_이름.R`의 함수 `slide_SNN()` + 문구 `text/ko/sNN.yaml`.
참고 구현: `slides/s09_scale.R`, `text/ko/s09.yaml`.

## 절대 규칙
1. **숫자를 손으로 쓰지 않는다.** 문구 파일(yaml)과 R 코드의 문자열에 결과 수치·설계값을 쓰지 않는다. 모든 수치는 아래 `d*` 함수나
   `f_*` 설계 상수 함수로 결과·config 파일에서 읽는다(출처가 추적표에 남는다). 허용되는 숫자는 식별자와 서지뿐이다
   (`text/ko/00_common.yaml`의 `check.number_contexts`: 저자 연도, NCT 번호, 표·그림 번호, "5~95백분위", "95% CI", "1:1", "2016 모델" 등).
   `check_deck.R`가 슬라이드 본문·표·발표자 노트의 모든 숫자를 그 슬라이드의 추적 행과 대조하고, 없으면 실패한다.
2. **근거는 시험 모집단(건강인, 체중 층화, B0)만.** 아토피 모집단·연구 범위 밖 체중(`results/atopic/`, `results/weight_generalization/`)은
   부록 A5에서만 쓴다. 본문 슬라이드에 "아토피", "환자 모집단"이라는 말도 쓰지 않는다.
3. **주분석 모형은 M1(체중 층 포함), M0 병기.** 경계 1종 오류는 `results/oc_models/type1_models.csv`(analysis_model 열)를 쓴다.
4. **제목은 그 슬라이드의 결론 문장**(26pt, 두 줄 이내). 위의 작은 kicker(13pt)는 논리 흐름 위치를 적는다(예: "논거 ② 판정 불안정").
5. **형식**: 요점(1단계 글머리표) 슬라이드당 5개 이하, 표 본문 6행 이하, 본문 16pt 이상, 표 12pt 이상. 넘치면 빌드가 멈춘다(추정 높이 포함).
6. **근거 태그**: `deck_slide(id, tag = "sim"|"lit"|"litsim"|"none")` → [모의], [문헌], [문헌+모의].
7. **발표자 노트**(`deck_notes`): 상세 설명, 주의 문구, 관련 보고서 절(예: "관련 보고서: 5.3절, Table 5-5"). 노트의 숫자도 추적 대상이다.
8. **약어**는 덱 전체에서 처음 나오는 슬라이드의 보이는 글자에 풀이를 둔다(`check.abbreviations`). 앞 슬라이드(S01~S08)가 대부분 정의한다.
   자기 슬라이드에서 처음 나올 수 있는 약어(TMDD, BLQ, ADA, IgG4, QC, BPD 등)는 풀어 쓴다: "표적 매개 약물 소실(TMDD)".
9. **대시 금지**: en-dash(–), em-dash(—), U+2212(−)를 쓰지 않는다. 범위는 수치 함수가 "a~b"로 만든다. 음수는 ASCII "-".
10. **사실만**: 해석이 필요한 문장은 결과로 뒷받침되는 것만 쓰고, 문장의 전제(예: "두 셀이 같은 시나리오")는 `premise(조건, "설명")`으로 검사한다.
11. 공용 파일(`lib/*`, `build_deck.R`, `check_deck.R`, `text/ko/00_common.yaml`, 다른 슬라이드 파일)은 고치지 않는다. 필요한 함수는
    자기 슬라이드 파일 안에 `sNN_` 접두어로 만든다. 공용 변경이 필요하면 결과 보고에 요청으로 적는다.

## 수치 함수 (lib/deck_lib.R; 모두 인쇄 문자열을 돌려주고 추적 행을 남김)
- `dv(rel, where, col, d, unit, item, scale)` 한 행의 값. `rel`은 `results/` 아래 상대 경로(예 "trialpop/tp_failure_by_set.csv") 또는 "config/..."
- `dint(rel, where, col, item, unit)` 정수(천 단위 쉼표). `dcount(rel, where, item)` 행 수.
- `drange(rel, where, col, d, unit, item, scale)` 여러 행의 최소~최대(두 모델 범위 등). 같으면 한 값.
- `dspan(rel, where, lo_col, hi_col, d, unit, item)` 여러 행의 lo 최솟값~hi 최댓값.
- `dext(rel, where, col, fun = max|min, d, unit, item)` 최댓값/최솟값.
- `dci(rel, where, est, lo, hi, d, unit, item)` "5.63% (95% CI 5.31~5.95)".
- `dcfg(file, path, item, fmt)` config 값. `dderived(item, rel, locator, raw, printed)` 파생값(차이·비·환산; locator에 계산식).
- `dsrc(item, rels)` 표·그림 전체 출처. `headline(p)` 표지 요약·논리 도식의 대표 수치(보고서/key numbers와 일치 검사 대상).
- 설계 상수(lib/deck_facts_common.R): `f_dose() f_n_arm() f_n_rand() f_wt_range() f_split() f_lloq() f_lloq_grid() f_limits() f_ci_level()
  f_nominal() f_study_days(schedule, "all"|"last"|"n") f_set(set, "r2"|"extrap"|"span") f_reps("boundary"|"ext")`.
- 행 조건 `where`는 data.table 식 문자열(그대로 추적표에 기록). 결과 파일 열 이름은 파일을 열어 확인한다. 보고서
  `regulatory/src/MS_report.Rmd`의 사실 청크가 같은 수치를 어떤 파일·조건으로 읽는지 보여 주므로 같은 조건을 쓰면 대조 검사가 쉬워진다.
- `item`은 영어로 짧게(추적표 설명). 한글 금지(추적표 영문 검사).

## 배치 함수 (단위 인치; 슬라이드 13.333 x 7.5; 여백 0.55)
- `deck_slide(id, tag, dark = FALSE)` 시작, `deck_end()` 끝(바닥글 자동: 출처 파일 ID, 버전, 커밋, 쪽 번호). 반드시 둘 다 부른다.
- `deck_kicker(s)`, `deck_title(s)`.
- `deck_bullets(items, box, size = 18)`: items 문자 벡터, "- "로 시작하면 2단계. `**굵게**`, `__강조색__` 표기 가능.
- `deck_text(s, box, size, bold, color, align, label, bg, geom)` 글상자(배경색·모양 가능: geom "roundRect").
- `deck_stat(value, label, box, color, bg)` 큰 수치 카드(40pt 수치 + 16pt 설명).
- `deck_box(box, fill, geom)` 도형(도식용; geom "roundRect", "rightArrow", "ellipse", "rect" 등 PowerPoint 도형 이름).
- `deck_table(df, box, widths, size = 13, highlight = 행번호)` 표(문자열 data.frame, 머리글 = names(df)).
- `deck_figure(p, name, box, src)` ggplot을 자리 크기 그대로 PNG로 그려 넣는다. `name`은 "sNN_설명"(겹치지 않게). `src`는 그림 자료 결과 파일.
  테마 `theme_deck(14)`, 색 `PAL$blue`, `PAL$orange`(두 모델: `MODEL_COL`, `MODEL_SHAPE`, `MODEL_LT`), 색 외에 표식·선 모양으로도 구분.
  축 글자 회전 금지(읽기 어렵다). 그림 안 글자도 문구 파일에서(`DK$txt$SNN$fig$...`). 그림 안 숫자(막대 값 표시)는 자료에서 바로 그린다.
  그림 제목은 넣지 않는다(슬라이드 제목이 결론). 보고서 PNG를 붙이지 않는다.
- 영역 상수: `GEO$ML`(0.55), `GEO$BODY_TOP`(1.65), `GEO$BODY_BOTTOM`(6.80), `GEO$CW`(12.23), `GEO$W`. 본문은 BODY_TOP~BODY_BOTTOM 안에.
- 배치는 슬라이드마다 다르게(좌 그림·우 요점, 카드 3개, 표+그림, 큰 수치 카드 등). 가로 전체 막대·가장자리 띠 장식 금지.

## 문구 파일 (text/ko/sNN.yaml)
최상위 키 = 슬라이드 ID(예 `S12:`). 흔한 키: `kicker`, `title`, `bullets`(목록), `table.head`, `fig.*`, `notes`(여러 줄 `|`).
자리표시자 `{이름}`은 `tx("S12.title", list(이름 = 값))`으로 채운다. 채워지지 않은 자리표시자는 빌드 오류.
문구는 한국어 본문 + 영문 기술 용어(AUC0-inf, λz, adjusted R², Cmax, GMR). 같은 키 구조로 `text/en/`을 만들면 영문판이 나온다.

## 개발 순환 (슬라이드 묶음마다 반복)
```
. <scratchpad>/renv_env.sh        # 또는 export RENV_CONFIG_EXTERNAL_LIBRARIES=... LANG=C.UTF-8 LC_ALL=C.UTF-8
Rscript reports/deck/build_deck.R --only S12,S13 --out <작업폴더>/g.pptx
Rscript reports/deck/check_deck.R <작업폴더>/g.pptx        # 5(약어) 실패는 부분 빌드에선 앞 슬라이드가 없어 생길 수 있다
DECK_SOFFICE_HELPER=<pptx skill>/scripts/office/soffice.py python3 reports/deck/render_deck.py <작업폴더>/g.pptx <작업폴더>/png --dpi 110
```
PNG를 직접 열어(이미지 보기) 글자 넘침·겹침·잘림, 표 가독성, 그림 해상도, 빈 값을 확인하고 고친다.
검사 9(렌더링 배치)는 렌더링한 PDF가 있어야 한다: render_deck.py를 먼저 실행한 뒤 check_deck.R를 실행한다(PDF는 pptx 옆이나 그 아래 폴더에서
같은 이름으로 찾고, 없으면 `--pdf 경로`). PDF의 글줄이 자기 도형 밖으로 2 pt 넘게 나가거나, 다른 도형의 글줄·그림과 겹치거나,
슬라이드 가장자리 0.3 in 안으로 들어오면 실패한다(`check_render.py`).
추정 높이: 줄 높이 = 1.2 x 줄 간격 배수 x 글자 크기(본문 줄 간격 1.1에서 실측 1.32배). 3% 넘게 넘치면 빌드가 멈춘다.
