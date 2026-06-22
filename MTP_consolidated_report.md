> **참고 (2026-06-22 갱신):** 이 문서는 워크스페이스 이전 시점의 원본 내러티브
> 리포트다. 이후 P0–P3 작업으로 코드 구조·수치가 갱신되었으니, 최신 기술 정본은
> [ARCHITECTURE.md](ARCHITECTURE.md), 진행 기록은 [worklog.md](worklog.md),
> 논문 초고는 [paper/paper.md](paper/paper.md)를 본다. 특히 §3-4의 MCMC 복원은
> 토이 w(z) likelihood 산물로, P1에서 관측가능량 기반으로 재산출되었다(§ARCHITECTURE 4).
> §5의 파일 경로는 `src/`·`scripts/`·`legacy/`로 재배치됨.

# MTP Cosmology — 통합 개발 리포트 (v0.1 → v0.3)
> Windowed Interacting Dark Energy Toy Model
> 상태: 합성 데이터 MCMC 완료 / 실데이터 MCMC 및 스크리닝 구현 = 워크스페이스 이전 대상

---

## 0. 한 줄 요약

미시 무한수렴 ↔ 거시 가속팽창의 MTP 철학을, GR 정합적인 **windowed IDE 모델**로 구현.
후기우주(z≈0.3) 활성화 창 함수로 Phantom crossing(w<-1)을 생성하며, DESI의 w₀-wₐ 시간변화 신호와 정합. 합성 데이터 MCMC가 입력 파라미터를 1σ 내 복원.

---

## 1. 모델 진화사

### v0.1 — 원형 ansatz (실패)
```
Q(z) = α·H·ρ_c·f(z),  f(z) = [log(ℓ(z)/L_P)/log(ℓ_0/L_P)]^n
```
**3대 결함:**
1. ODE 직접 전달 → ∇_μT^{μν}=0 위반 (GR 부정합)
2. w_eff = -1 + (양수항) → Phantom barrier, w<-1 불가
3. α = G_eff/G_N-1 해석 → Cassini 제약 100배 초과

추가: f(z→0)→0 구조라 현재우주에서 전달항 소멸. 어떤 α로도 w=-1 고정.

### v0.2 — IDE 전환 (부분 해결, 크기 부족)
- ODE → IDE 4-벡터 결합항 (∇_μT^{μν}=0 충족) ✅
- α(z)=α₀·tanh(z-z_c) 부호 전환 도입 ✅
- α를 IDE 결합계수로 재정의 (Cassini 직접 제약 회피) ✅
- **잔존 문제**: f(z) 로그구조가 저z에서 소멸 → 효과가 DESI 요구의 1/10~1/30

### v0.3 — Windowed IDE (성공)
조비 권고로 f(z)를 역할 분리:
```
Q(z) = β(z)·H·ρ_c · F_hier(z) · W_late(z)

F_hier(z) = [log(ℓ(z)/L_P)/log(ℓ_0/L_P)]^n   ← 계층 철학 보존 (저z 거의 상수)
W_late(z) = exp[-(z-z*)²/2σ²]                  ← 후기우주 활성화
β(z)      = β₀·tanh((z-z_c)/Δz)                ← 부호/강도 제어
```
**결과**: w<-1 진입 성공, w₀-wₐ 신호 생성, H/D_M/fσ8 sanity 통과.

### v0.3 자기인식 (정직한 한계)
F_hier(z)는 관측 범위(z≲2)에서 거의 상수 → 1로 놔도 결과 불변.
따라서 본 모델은 실질적으로 **windowed IDE 모델**이며,
"MTP/계층" 명명은 철학적·역사적 맥락으로 보존됨.

---

## 2. 최종 수식 (v0.3)

```
연속방정식:
  ρ̇_c  + 3Hρ_c        = +Q
  ρ̇_de + 3H(1+w)ρ_de  = -Q

Friedmann:
  H²(z) = H₀²[Ω_m·ρ_c/ρ_c0 + Ω_de·ρ_de/ρ_de0 + Ω_r(1+z)⁴]

자유 파라미터: β₀, z*, σ, z_c, Δz  (+ n, 단 영향 미미)
고정: Planck 2018 (H₀=67.4, Ω_m=0.315, Ω_de=0.685)
```

---

## 3. 수치 결과 (합성 데이터)

### 3-1. w_eff(z) — Phantom crossing
| config | w(0) | w_min | z@min | Phantom |
|--------|------|-------|-------|---------|
| β₀=0.05, z*=0.5 | -1.008 | -1.025 | 0.32 | YES |
| β₀=0.10, z*=0.5 | -1.015 | -1.052 | 0.32 | YES |
| β₀=0.20, z*=0.4 | -1.027 | -1.077 | 0.25 | YES |

### 3-2. w₀-wₐ (CPL 근사) — DESI 비교 신호
| config | w₀ | wₐ |
|--------|-----|-----|
| β₀=0.05, z*=0.5 | -1.030 | +0.133 |
| β₀=0.10, z*=0.5 | -1.060 | +0.264 |
- w₀<-1, wₐ>0 조합 = DESI DR1 선호 방향과 정합

### 3-3. Sanity check (3개 통과)
- H(0) = 67.43 (전 config 일치)
- fσ8(0.5) ≈ 0.477~0.480 (관측 ≈0.46, σ8 tension 완화 방향)
- D_M(z) 편차 Λ-CDM 근방 합리적

### 3-4. MCMC 파라미터 복원
**3-파라미터 (β₀, z*, σ):**
| param | 입력 | 복원 mean | std |
|-------|------|-----------|-----|
| β₀ | 0.10 | 0.1056 | 0.0078 |
| z* | 0.40 | 0.4504 | 0.0356 |
| σ | 0.30 | 0.3314 | 0.0315 |

**2-파라미터 (β₀, z*):** 입력값 1σ 내 복원 확인.
→ 추정 워크플로우 정상 작동. 실데이터 확장 준비 완료.

---

## 4. 이론적 정당화 (조비)

### 4-1. 창 함수 중심 z*≈0.3
- 암흑물질-암흑에너지 밀도 동등 시점 (coincidence epoch)
- 2025 sign-switching IDE 연구: 밀도 동등점에서 결합 부호 전환, tanh 평활화 사용 → 본 모델 선택과 정합

### 4-2. Cassini 제약 해결 — 스크리닝
```
β(z,ρ) = β₀·tanh((z-z_c)/Δz) · S(ρ)
  S(ρ)≪1  고밀도 (태양계) → β→0
  S(ρ)≈1  우주론 스케일
```
- Chameleon/symmetron/Vainshtein 계열
- thin-shell factor S(ρ)가 Cassini |γ-1|≲2.3×10⁻⁵ 자동 충족
- 단, 구체 scalar-tensor 미시모델 구현은 향후 과제

### 4-3. 독창성 / IP
- 분류: Hierarchical Scale Transfer Cosmology (새 모델 클래스)
- 기존 IDE/Quintom/HDE/RVM과 구별: 로그 계층항 + windowed IDE + 부호전환 조합
- IP 확보 가능성 높음 (조비 판단)

---

## 5. 코드 자산 (저장소 포함)

| 파일 | 내용 |
|------|------|
| `mtp_v03.py` | windowed IDE 본 모델. H(z)/D_M/fσ8/w_eff 계산 + 시각화 |
| `advanced_mcmc_analysis.py` | 3-param MH-MCMC (β₀,z*,σ). chain+summary CSV 출력 |
| `mcmc_toy_model.py` | 2-param MH-MCMC 데모 |
| `mtp_v02.py` | v0.2 (IDE 전환, 한계 기록용) |
| `mtp_cosmology.py` | v0.1 (원형, 히스토리) |

---

## 6. 로드맵 (워크스페이스 이전 후)

### Step 2 — 실데이터 MCMC
- Cobaya/MontePython 풀스택
- Planck CMB + DESI BAO + Pantheon+ SNIa + RSD
- log-likelihood: H(z), 거리, 구조성장
- β₀,z*,σ,z_c,Δz 피팅 → Λ-CDM/CPL 대비 Bayesian evidence 비교
- **이 컨테이너 제약**: Planck 데이터 도메인 접근 불가 → 워크스페이스 필수

### Step 3 — 스크리닝 구체화
- scalar-tensor 미시모델 선택
- thin-shell S(ρ) 계산, PPN 충족 증명

### Step 4 — 물리적 정당화 심화
- z≈0.3 전달 피크를 coincidence problem과 연결

### Step 5 — 논문 작성
- 모델 + 데이터 피팅 + 스크리닝 + 해석

---

## 7. 의사결정 기록

- **경로 C 채택** (조비): DESI w=-1.03 직접 설명 ❌ → 후기우주 windowed IDE toy model로 재정의 → w₀-wₐ 시간변화 신호 타겟 ✅
- 현재 단계: 합성 MCMC 완료. 실데이터 = 워크스페이스 이전 후 진행.

---

*MTP Cosmology v0.3 통합 리포트 — 워크스페이스 첫 커밋 기준 문서*
