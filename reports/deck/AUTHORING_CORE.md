# 핵심 덱 v1.1 작성 규칙 (reports/deck, `--deck core`)

핵심 덱은 결과보고 덱(v1.0.1, 32장, 기술 백업)과 같은 도구(R officer, `lib/deck_lib.R`)로 만든다. 공통 규칙은 `AUTHORING.md`를 따르고,
아래는 핵심 덱에서 다르거나 더해진 것이다(지시 2026-09-29, 사전 등록 `config/prereg_20260929.yaml`, DECISIONS D-063).

## 파일
- 슬라이드 코드: `slides_core/sN_*.R`(함수 `slide_SN()`), 별첨 `slides_core/aN_*.R`(함수 `slide_AN()`, 예 `slide_A3a()`).
- 문구: `text/ko_core/sN.yaml`, `aN.yaml`(최상위 키 = 슬라이드 ID). 순서·약어·숫자 허용 문맥: `text/ko_core/00_common.yaml`.
- 공용 도우미: `slides_core/c00_helpers.R`(`core_title`, `core_body`, `core_caption`, `core_card`, `core_shaded_panels`, `core_sap`, `fl`/`cl` 내림·올림).
- 빌드: `Rscript reports/deck/build_deck.R --deck core [--only S5,A3a --out <경로>.pptx]`. 검사: `Rscript reports/deck/check_deck.R <pptx> --deck core`(렌더링 뒤).
- 렌더링: `python3 reports/deck/render_deck.py <pptx> <PNG 폴더> --dpi 100`(PDF도 같은 폴더에).

## 형식(빌드와 검사 6이 강제)
- 글자: 제목 28pt 굵게(2줄 이내, 넘으면 빌드 중단), 킥커 16pt, 본문 18pt, 표·그림 안 글자·캡션·각주 14pt 이상, 바닥글 10pt.
- 그림 글자: `theme_core(16)`을 쓰고 `geom_text`/`annotate` 크기는 `PT(14)` 이상(`PT(pt)`는 pt를 ggplot 크기로 바꾼다). `deck_figure`가 테마와 글자 레이어를 검사해 14pt 미만이면 멈춘다.
- 본문 슬라이드(S2~S10): 본문(이름이 `body`로 시작하는 글상자) 3줄 이내, 주 시각 요소(그림·도식·카드 묶음) 면적이 본문 영역의 60% 이상.
- 요점(글머리표) 슬라이드당 3개 이하, 표 본문 8행 이하. 별첨은 60%·3줄 한도가 없지만 한 장 한 메시지, 글자 하한은 같다.
- 색: AUClast 파랑(`PAL$blue`), 외삽·AUCinf 주황(`PAL$orange`). 모델 구분은 색이 아니라 표식·선 모양(`CORE_MODEL_SHAPE`, `CORE_MODEL_LT`).
- 바닥글: 출처 파일 ID와 버전만(자동). 근거 태그 `deck_slide(id, tag = "sim"|"lit"|"litsim"|"none")`.

## 용어(고정)
- 신뢰할 수 있는 AUCinf = adjusted R² ≥ 0.90 그리고 외삽 ≤ 20%(기준 세트 (iii)). 수치는 `f_set("iii", "r2")`, `f_set("iii", "extrap")`.
- AUClast/AUCinf = 참값 기준(관측 tlast까지의 참 AUC / 모델 적분 참 AUCinf; 결과 열 coverage_true 또는 "window coverage").
- 2016 모델(주 모델), 2020 모델(민감도 모델). 짧게는 `DK$txt$common$models_short`.
- 본문(S1~S10)에서 쓰지 않는 말: Wilson 분류, M0/M2, 세트 (ii)/(iv), 규칙 B/C, 칸 수 세부, 모델 파라미터. 별첨에서는 쓸 수 있다(풀어서).
- 결과 파일의 옛 표기(AUC0-last, AUC0-inf)는 덱 글자에서 AUClast, AUCinf로 쓴다. 다른 슬라이드는 "슬라이드 N"이 아니라 ID(S5, 별첨 A3a)로 가리킨다.

## 수치와 약어
- 숫자는 문구 파일에 쓰지 않는다(`{자리표시자}`). 모든 수치는 `d*` 함수·`f_*` 설계 상수로 결과·config 파일에서 읽는다. 발표자 노트의 숫자도 검사 대상이다.
- 약어는 덱 전체에서 처음 보이는 슬라이드의 글자에 풀이가 있어야 한다(표지는 용어 줄로 대신). S1~S10에서 이미 풀이한 약어: PK, AUClast, AUCinf, Cmax, NCA, λz, LLOQ, FDA, BLA, SAP, NCT, tlast, EMA.
  별첨은 순서가 바뀔 수 있으므로 그 밖의 약어(GMR, CI, CV, Km, Vmax, V2, ka, ke, ADA, BLQ, SD, SE, TMDD, QC, ICH, Clast, tmax 등)는 각 별첨 슬라이드에서 처음 쓸 때 풀어 쓴다(예: "기하평균비(GMR)").
- 아토피·체중 범위 밖 모집단 결과(`results/atopic/`, `results/weight_generalization/`)는 별첨 A11에서만 쓴다.
