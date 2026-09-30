# 핵심 덱 작성 규칙 (v1.2; reports/deck, `--deck core`)

핵심 덱은 결과보고 덱(v1.0.1, 32장, 기술 백업)과 같은 도구(R officer, `lib/deck_lib.R`)로 만든다. 공통 규칙은 `AUTHORING.md`를 따르고,
아래는 핵심 덱에서 다르거나 더해진 것이다(지시 2026-09-29, 사전 등록 `config/prereg_20260929.yaml`·`config/prereg_20260929_oc.yaml`, DECISIONS D-063~D-065).

## v1.2 구조(D-065)
- 본문 8장: S1 표지, S2 결론, S3 방법, S4 전제, S5 ①, S6 ②, S7 ③(1종·2종 오류), S8 결론·제안. 본문 한도(`main_limits`)는 S2~S8.
- 별첨 A1~A14(A13 없음), 보강 B1(표시 "보강 A-①"), B2("보강 A-②"). ③의 세부는 A5a~A5f. 순서는 `00_common.yaml`의 `order`.
- 본문 슬라이드는 관련 별첨 한 줄을 문구 파일의 `xref`에 둔다(`deck_slide`가 본문 영역 아래 14pt로 그리고 본문 아래 끝을 올린다). S2 카드와 S8 근거 줄에는 본문 수치만 쓴다.
- S7·S2 카드 3·S8 근거 ③의 문구는 사전 등록 문구 규칙(`results/oc_curves/oc_wording_rules.csv`)으로 고른다. 규칙이 달라지면 문구가 따라 바뀌므로 문구 파일에 대안 문장을 함께 둔다.

## 작업 방식(초안 모드, 사용자 지시 2026-09-29, D-065)
- 덱 수정은 기본으로 초안 모드: 바뀐 슬라이드만 `--only`로 빌드·검사(숫자 추적과 검사 9종 유지), 전체 렌더는 마지막에 1회.
- 독립 검토 에이전트는 사용자가 요청할 때만 쓴다. 규제 패키지(scripts/60)는 사용자가 "동결"을 지시할 때만 다시 만든다. 그 전까지 check_deck 8의 자기 이전 게시 행과의 불일치는 PENDING FREEZE(통과)로 표시된다.
- 새 모의가 필요하면 시작 전에 예상 소요 시간을 보고하고, 1시간을 넘으면 조건당 약 500회 예비 실행으로 방향을 먼저 확인한다. 장시간 실행은 배치 단위로 이어서 돌릴 수 있게 한다(scripts/64 방식).

## 파일
- 슬라이드 코드: `slides_core/sN_*.R`(함수 `slide_SN()`), 별첨 `slides_core/aN_*.R`(함수 `slide_AN()`, 예 `slide_A3a()`).
- 문구: `text/ko_core/sN.yaml`, `aN.yaml`(최상위 키 = 슬라이드 ID). 순서·약어·숫자 허용 문맥: `text/ko_core/00_common.yaml`.
- 공용 도우미: `slides_core/c00_helpers.R`(`core_title`, `core_body`, `core_caption`, `core_card`, `core_shaded_panels`, `core_sap`, `fl`/`cl` 내림·올림).
- 빌드: `Rscript reports/deck/build_deck.R --deck core [--only S5,A3a --out <경로>.pptx]`. 검사: `Rscript reports/deck/check_deck.R <pptx> --deck core`(렌더링 뒤).
- 렌더링: `python3 reports/deck/render_deck.py <pptx> <PNG 폴더> --dpi 100`(PDF도 같은 폴더에).

## 형식(빌드와 검사 6이 강제)
- 글자: 제목 28pt 굵게(2줄 이내, 넘으면 빌드 중단), 킥커 16pt, 본문 18pt, 표·그림 안 글자·캡션·각주 14pt 이상, 바닥글 10pt.
- 그림 글자: `theme_core(16)`을 쓰고 `geom_text`/`annotate` 크기는 `PT(14)` 이상(`PT(pt)`는 pt를 ggplot 크기로 바꾼다). `deck_figure`가 테마와 글자 레이어를 검사해 14pt 미만이면 멈춘다.
- 본문 슬라이드(S2~S8): 본문(이름이 `body`로 시작하는 글상자) 3줄 이내, 주 시각 요소(그림·도식·카드 묶음) 면적이 본문 영역의 60% 이상.
- 요점(글머리표) 슬라이드당 3개 이하, 표 본문 8행 이하. 별첨은 60%·3줄 한도가 없지만 한 장 한 메시지, 글자 하한은 같다.
- 표: `deck_table`(8행 이하). 셀 위아래 여백은 `pad`(기본 4pt, 자리가 모자라면 2pt), 묶음 구분은 `group_end`(그 행들 아래에만 굵은 선). 상자 높이는 `deck_table_h`로 맞춘다. LibreOffice는 열 폭에 가까운 셀의 행을 높게 그리므로 셀 글자는 짧게 둔다(숫자와 한글 사이의 줄바꿈 방지 문자도 행을 높일 수 있다).
- 색: AUClast 파랑(`PAL$blue`), 외삽·AUCinf 주황(`PAL$orange`). 모델 구분은 색이 아니라 표식·선 모양(`CORE_MODEL_SHAPE`, `CORE_MODEL_LT`).
- 바닥글: 출처 파일 ID와 버전만(자동). 근거 태그 `deck_slide(id, tag = "sim"|"lit"|"litsim"|"none")`.

## 용어(고정)
- 신뢰할 수 있는 AUCinf = adjusted R² ≥ 0.90 그리고 외삽 ≤ 20%(기준 세트 (iii)). 수치는 `f_set("iii", "r2")`, `f_set("iii", "extrap")`.
- AUClast/AUCinf = 참값 기준(관측 tlast까지의 참 AUC / 모델 적분 참 AUCinf; 결과 열 coverage_true 또는 "window coverage").
- 2020 모델(주 모델), 2016 모델(민감도 모델)(v1.2, 사용자 지시, D-065). 그림·표·문장 모두 2020 먼저. 짧게는 `DK$txt$common$models_short`(순서 k2020, k2016). 모델 표식 `CORE_MODEL_SHAPE`(2020 검은 원, 2016 회색 삼각형).
- 판정 구성: AUClast + Cmax(제안), + AUCinf(세 지표), AUCinf + Cmax(가이드라인 기본). 별첨에서 규칙을 밝힐 때는 AUCinf(A) = 신뢰할 수 있는 AUCinf만, AUCinf(B) = λz 산출 전원. 1종 오류 = 참 AUCinf 비가 동등 한계인 경계 칸에서 동등 판정 비율, 2종 오류 = 동등한 제품(참 AUCinf 비와 참 Cmax 비가 모두 한계 안)을 동등으로 판정하지 못할 확률. "1 쪽 치우침" = 추정 기하평균비가 참 비보다 1에 가까운 쪽으로 간 양.
- 본문(S1~S8)에서 쓰지 않는 말: Wilson 분류, M0/M2, 세트 (ii)/(iv), 규칙 B/C, 칸 수 세부, 모델 파라미터 약어(Vmax 등). 별첨에서는 쓸 수 있다(풀어서). 예외: S7 본문의 "신뢰 기준 적용 시/적용하지 않아도"는 규칙 A/B를 풀어 쓴 말이다.
- 대상자 구분은 "기준 미달자 / 충족자"(신뢰할 수 있는 AUCinf 기준), 두 투여군은 "arm"(예: "arm당", "두 arm 모두"), 모델 값과 비교하는 값은 "참값"(실제 채혈 시각·실제 Phoenix처럼 관측·실물을 가리킬 때만 "실제").
- 1종 오류 칸 분류: "초과"는 Wilson 분류(95% 구간 하한 > 5%)에만 쓰고, 점추정 비교는 "점추정 5% 초과"로 적는다(본문 S2~S8은 점추정만, Wilson 분류는 별첨 A5a·A5b). M1은 "주분석 후보"(의뢰자 결정 전).
- 모의 최솟값은 추출마다 변한다: "최솟값"을 적을 때는 한 번의 모의(모델별 20,000명) 값임을 문구나 노트에 밝힌다(D-063).
- 결과 파일의 옛 표기(AUC0-last, AUC0-inf)는 덱 글자에서 AUClast, AUCinf로 쓴다. 다른 슬라이드는 "슬라이드 N"이 아니라 ID(S5, 별첨 A3a)로 가리킨다.

## 수치와 약어
- 숫자는 문구 파일에 쓰지 않는다(`{자리표시자}`). 허용 문맥(`check.number_contexts`)에 "1 쪽", "규칙 4a~4c", 별첨 ID(A5a~A5f 등)가 있다. 모든 수치는 `d*` 함수·`f_*` 설계 상수로 결과·config 파일에서 읽는다. 발표자 노트의 숫자도 검사 대상이다.
- 약어는 덱 전체에서 처음 보이는 슬라이드의 글자에 풀이가 있어야 한다(표지는 용어 줄로 대신). S1~S8에서 이미 풀이한 약어: PK, AUClast, AUCinf, Cmax, NCA, λz, LLOQ, FDA, BLA, SAP, NCT, tlast, EMA.
  별첨은 순서가 바뀔 수 있으므로 그 밖의 약어(GMR, CI, CV, Km, Vmax, V2, ka, ke, ADA, BLQ, SD, SE, TMDD, QC, ICH, Clast, tmax 등)는 각 별첨 슬라이드에서 처음 쓸 때 풀어 쓴다(예: "기하평균비(GMR)").
- 아토피·체중 범위 밖 모집단 결과(`results/atopic/`, `results/weight_generalization/`)는 별첨 A11에서만 쓴다.
