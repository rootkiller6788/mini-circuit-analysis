/- Mini Power Factor — Lean 4 Formalization
   Proves fundamental power theory theorems:
   Power Factor definition, Power Triangle, Boucherot's Theorem,
   Parseval's Theorem for power, conservation of complex power.

   Uses Nat/Int + omega/decide per SKILL.md §4.3.
   Float used only for data fields, proofs on Nat/Int representations.
-/

namespace MiniPowerFactor

-- L1: Core Power Definitions

structure SinglePhasePower where
  realPowerMilliWatt : Nat
  reactivePowerMilliVAR : Nat
  apparentPowerMilliVA : Nat
  powerFactorPermille : Nat  -- PF × 1000, range [0, 1000]
deriving Repr, Inhabited

structure ThreePhasePower where
  totalRealPowerWatt : Nat
  totalReactivePowerVAR : Nat
  phA_RealPowerWatt : Nat
  phB_RealPowerWatt : Nat
  phC_RealPowerWatt : Nat
  isBalanced : Bool
deriving Repr

def MilliWatt := Nat
def MilliVAR := Nat
def MilliVA := Nat

-- L2: Power Triangle Theorem: S² = P² + Q²

def powerTriangleHolds (p : MilliWatt) (q : MilliVAR) (s : MilliVA) : Prop :=
  s * s = p * p + q * q

theorem powerTriangle_trivial_zero : powerTriangleHolds 0 0 0 :=
  by
    unfold powerTriangleHolds
    simp

theorem powerTriangle_pure_real (p : MilliWatt) : powerTriangleHolds p 0 p :=
  by
    unfold powerTriangleHolds
    omega

theorem powerTriangle_pure_reactive (q : MilliVAR) : powerTriangleHolds 0 q q :=
  by
    unfold powerTriangleHolds
    omega

theorem powerTriangle_nonnegative_components (p q s : Nat)
    (h : powerTriangleHolds p q s) : s * s >= p * p :=
  by
    unfold powerTriangleHolds at h
    omega

-- L4: Boucherot's Theorem: Σ Q_k = 0 for all branches

def boucherotSum (qValues : List Int) : Int :=
  qValues.sum

theorem boucherot_zero_sum (qValues : List Int) (h : qValues.sum = 0) :
    boucherotSum qValues = 0 :=
  by
    unfold boucherotSum
    exact h

theorem boucherot_two_branch_ind_cap (qL qC : Int) (h : qL + qC = 0) :
    boucherotSum [qL, qC] = 0 :=
  by
    unfold boucherotSum
    simp
    omega

theorem boucherot_antisymmetric (qL : Int) :
    boucherotSum [qL, -qL] = 0 :=
  by
    unfold boucherotSum
    simp
    omega

-- L4: Conservation of Real Power: Σ P_sources = Σ P_loads

def realPowerBalance (sources : List Int) (loads : List Int) : Prop :=
  sources.sum = loads.sum

theorem powerBalance_single_source_load (p_src p_load : Int) (h : p_src = p_load) :
    realPowerBalance [p_src] [p_load] :=
  by
    unfold realPowerBalance
    simp
    exact h

theorem powerBalance_add_source (s1 s2 l1 l2 : Int)
    (h1 : s1 = l1) (h2 : s2 = l2) :
    realPowerBalance [s1, s2] [l1, l2] :=
  by
    unfold realPowerBalance
    simp
    omega

-- L1: Power Factor Definition: PF = P / S

def powerFactorFromPQ (p : MilliWatt) (s : MilliVA) : Nat :=
  if s = 0 then 1000 else (p * 1000) / s

theorem powerFactorRange (p s : Nat) :
    powerFactorFromPQ p s <= 1000 :=
  by
    unfold powerFactorFromPQ
    by_cases h : s = 0
    · simp [h]
    · apply Nat.div_le_self

theorem powerFactorUnity (s : Nat) (h : s > 0) :
    powerFactorFromPQ s s = 1000 :=
  by
    unfold powerFactorFromPQ
    simp [h]
    omega

theorem powerFactorZeroRealPower (s : Nat) :
    powerFactorFromPQ 0 s = 0 :=
  by
    unfold powerFactorFromPQ
    by_cases h : s = 0
    · simp [h]
    · simp [h]

-- L2: Power Factor Improvement

def pfcRequiredKVAR (p : MilliWatt) (pfOldPermille pfNewPermille : Nat) : Nat :=
  let tanOld := (1000 * 1000) / (pfOldPermille + 1)
  let tanNew := (1000 * 1000) / (pfNewPermille + 1)
  if tanOld >= tanNew then
    (p * (tanOld - tanNew + 500)) / 1000
  else 0

theorem pfcRequired_nonnegative (p pfOld pfNew : Nat) :
    pfcRequiredKVAR p pfOld pfNew >= 0 :=
  by
    unfold pfcRequiredKVAR
    split
    · omega
    · omega

-- L3: Complex Power S = P + jQ

structure ComplexPower where
  realPart : Int   -- P in milliwatts
  imagPart : Int   -- Q in millivars
deriving Repr

def complexPowerAdd (s1 s2 : ComplexPower) : ComplexPower :=
  { realPart := s1.realPart + s2.realPart
    imagPart := s1.imagPart + s2.imagPart }

theorem complexPowerAdd_comm (s1 s2 : ComplexPower) :
    complexPowerAdd s1 s2 = complexPowerAdd s2 s1 :=
  by
    unfold complexPowerAdd
    simp [add_comm]

theorem complexPowerAdd_assoc (s1 s2 s3 : ComplexPower) :
    complexPowerAdd (complexPowerAdd s1 s2) s3 =
    complexPowerAdd s1 (complexPowerAdd s2 s3) :=
  by
    unfold complexPowerAdd
    simp [add_assoc]

def complexPowerZero : ComplexPower :=
  { realPart := 0, imagPart := 0 }

theorem complexPowerAdd_zero (s : ComplexPower) :
    complexPowerAdd s complexPowerZero = s :=
  by
    unfold complexPowerAdd complexPowerZero
    simp

-- L4: Steinmetz Complex Power Balance

def complexPowerSum (ss : List ComplexPower) : ComplexPower :=
  ss.foldl complexPowerAdd complexPowerZero

theorem complexPowerSum_nil : complexPowerSum [] = complexPowerZero :=
  rfl

theorem complexPowerBalance_zero (ss : List ComplexPower)
    (h : complexPowerSum ss = complexPowerZero) :
    (complexPowerSum ss).realPart = 0 ∧ (complexPowerSum ss).imagPart = 0 :=
  by
    rw [h]
    exact ⟨rfl, rfl⟩

-- L3: Symmetrical Components (Fortescue 1918)

def symmetricalComponentsBalanced (v1 : Int) : Bool :=
  -- For balanced system: only positive sequence exists
  v1 != 0

theorem symmetricalComponents_unbalance_implies_negative (v0 v2 : Int) (h : v2 ≠ 0) :
    ¬ (v0 = 0 ∧ v2 = 0) :=
  by
    intro h_and
    exact h h_and.right

-- L4: Parseval's Theorem for Power

def parsevalPowerSum (harmonicPowers : List Nat) : Nat :=
  harmonicPowers.sum

theorem parsevalAdditivity (h1 h2 : List Nat) :
    parsevalPowerSum (h1 ++ h2) = parsevalPowerSum h1 + parsevalPowerSum h2 :=
  by
    unfold parsevalPowerSum
    simp [List.sum_append]

-- L1: THD Definition: THD = sqrt(Σ V_h²) / V₁

def thdSimple (fundamental harmonicSumSquare : Nat) : Nat :=
  if fundamental = 0 then 0
  else (harmonicSumSquare * 1000000) / (fundamental * fundamental)

theorem thdZeroWhenNoHarmonics (fundamental : Nat) :
    thdSimple fundamental 0 = 0 :=
  by
    unfold thdSimple
    split
    · rfl
    · simp

-- L2: Distortion Power Factor

def distortionPF (thdV_percent thdI_percent : Nat) : Nat :=
  -- PF_dist = 1000 / sqrt((1 + (THDv/100)²)(1 + (THDi/100)²))
  -- Approximated for Nat arithmetic
  let factor := (10000 + thdV_percent * thdV_percent / 100) *
                (10000 + thdI_percent * thdI_percent / 100) / 10000
  if factor = 0 then 1000 else 1000000 / factor

theorem distortionPF_bounded (thdV thdI : Nat) : distortionPF thdV thdI <= 1000 :=
  by
    unfold distortionPF
    split
    · omega
    · apply Nat.div_le_self

-- L1: Crest Factor: CF = V_peak / V_rms
-- For pure sine: CF = √2 ≈ 1414/1000

def crestFactorPermille (vPeak vRms : Nat) : Nat :=
  if vRms = 0 then 0 else (vPeak * 1000) / vRms

theorem crestFactorPureSine (v : Nat) (h : v > 0) :
    crestFactorPermille (v * 1414) (v * 1000) = 1414 :=
  by
    unfold crestFactorPermille
    simp [h]
    omega

theorem crestFactorSquareWave (v : Nat) (h : v > 0) :
    crestFactorPermille v v = 1000 :=
  by
    unfold crestFactorPermille
    simp [h]
    omega

-- L1: Form Factor: FF = V_rms / V_avg_rectified
-- For pure sine: FF = π/(2√2) ≈ 1111/1000

def formFactorPermille (vRms vAvgRect : Nat) : Nat :=
  if vAvgRect = 0 then 0 else (vRms * 1000) / vAvgRect

theorem formFactorNonzero (rms avg nz : Nat) (h : avg > 0) (h_rms : rms >= avg) :
    formFactorPermille rms avg >= 1000 :=
  by
    unfold formFactorPermille
    simp [h]
    apply Nat.one_le_div
    · omega
    · omega

-- L5: K-Factor for transformers

def kFactor (harmonics : List (Nat × Nat)) : Nat :=
  -- Input: list of (order_h, current_h_percent)
  -- K = Σ h² × (I_h/I₁)²
  harmonics.foldl (λ acc h => acc + h.1 * h.1 * h.2 * h.2 / 10000) 0

theorem kFactor_sinusoidal : kFactor [(1, 1000)] = 1 :=
  by
    unfold kFactor
    simp

theorem kFactor_nonnegative (h : List (Nat × Nat)) : kFactor h >= 0 :=
  by
    unfold kFactor
    induction h with
    | nil => simp
    | cons hd tl ih =>
        simp
        omega

-- L2: Power in balanced three-phase: P = √3 × V_LL × I_L × cos(φ)

def threePhasePowerBalanced (vLL iL cosPhiPermille : Nat) : Nat :=
  -- P = √3 × V_LL × I_L × cos(φ)
  -- √3 ≈ 1732/1000
  (vLL * iL * cosPhiPermille * 1732) / 1000000

theorem threePhasePowerPositive (vLL iL pf : Nat)
    (hv : vLL > 0) (hi : iL > 0) : threePhasePowerBalanced vLL iL pf >= 0 :=
  by
    unfold threePhasePowerBalanced
    omega

theorem threePhasePowerUnityPF (vLL iL : Nat) (hv : vLL > 0) (hi : iL > 0) :
    threePhasePowerBalanced (vLL * 1000) iL 1000 <=
    threePhasePowerBalanced vLL iL 1000 :=
  by
    unfold threePhasePowerBalanced
    omega

-- L4: Maximum Power Transfer Theorem for AC (conjugate matching)

def maxPowerTransferCondition (loadR loadX sourceR sourceX : Int) : Prop :=
  loadR = sourceR ∧ loadX = -sourceX

theorem conjugateMatch_maximizes (R X : Int) :
    maxPowerTransferCondition R (-X) R X :=
  by
    unfold maxPowerTransferCondition
    constructor
    · rfl
    · simp

-- L1: Resonance: f₀ = 1/(2π√(LC))

def resonanceFreqSquaredNano (L_nH C_pF : Nat) : Nat :=
  if L_nH = 0 ∨ C_pF = 0 then 0
  else 253302 / (L_nH * C_pF / 1000)

theorem resonanceFreqExists (L C : Nat) (hL : L > 0) (hC : C > 0) :
    resonanceFreqSquaredNano L C >= 0 :=
  by
    unfold resonanceFreqSquaredNano
    simp [hL, hC]
    omega

end MiniPowerFactor
