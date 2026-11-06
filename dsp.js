// SPDX-License-Identifier: CC0-1.0
/*
 * x0xb0x VCF in real time: a model of the ngspice model in ../vcf-simulation.
 *
 * The same file is the AudioWorklet that index.html plays (registerProcessor at the end)
 * and the engine of the offline renders (render.mjs, through module.exports). Every value
 * is a part value from the KiCad schematic (the vcf-parts.json the ngspice deck reads), a
 * device number from the study's models, a number the study states, or one of the fitted
 * constants in FIT, each commented with what it was fitted to. See README.md.
 *
 * Signal path, all voltages deviations from the operating point:
 *   VCO (N$126) -> R62/C17/R70 -> Q12 differential pair (tanh) -> four-rung diode ladder
 *   (C18 at the bottom, C19, C24, C26; each diode pair I*tanh(dv/2VT)) -> C25/C27 -> Q21
 *   differential amplifier (tanh) -> Q19 -> C14 -> Resonance pot VR4a -> C15 -> Q18 -> C23
 *   -> R97/R98 -> back to Q12's other base. The output tap C20/R122, C21/R121 gives the
 *   VCA input N$623. The ladder is solved implicitly (trapezoidal rule, Newton on four
 *   unknowns with the whole feedback loop inside), so there is no unit delay in the loop.
 *
 * Control path, at the host rate, in absolute volts:
 *   Cutoff VR3 -> TM3 -> R73 sets a current; Q10/Q11 make the ladder current
 *   I = K * (I(R73) - Q11's base current) * exp((V134 - V561)/VT), where N$561 is Q9's
 *   emitter. The summing node N$134 gets Env Mod (MEG through VR5 and R63) and Accent (MEG
 *   through the 4066, VR7, D24, R46, VR4b with C13, R72). The MEG is C62, struck through D37
 *   and discharged through R138 + VR6 into a diode stack.
 */
'use strict';

const VT = 0.025693;          // kT/q at 25 C, the temperature of the ngspice runs
const IS_D = 2.52e-9;         // 1N4148 as in the study's deck: IS 2.52 nA, N 1.752
const NVT_D = 1.752 * VT;
const INV2VT = 1 / (2 * VT);

// Part values (vcf-parts.json, fullchain-parts.json) and the study's trims.
const PART = {
  // Ladder and its feedback amplifier
  C: 33e-9, C18: 18e-9, C25: 0.1e-6, R114: 100e3, R110: 22e3,
  C14: 1e-6, C15: 1e-6, VR4: 50e3, R67: 2.2e3, R99: 100e3, R66: 100e3, C22: 1e-6, R94: 10e3, C23: 1e-6, R97: 10e3, R98: 2.2e3, R96: 10e3,
  C20: 0.01e-6, C21: 0.01e-6, R121: 220e3, R122: 100e3, R124: 2.2e3,
  R62: 220e3, R70: 2.2e3, C17: 1e-6, R68: 2.2e3, C16: 10e-6,
  // Control
  // N$561 at rest (DC table), and Q9's base divider R64 || R65
  V561: 3.26, R64_R65: 5e3, R71: 2.2e3, R63: 220e3, R61: 10e3, VR5: 50e3, R72: 100e3, R46: 47e3, C13: 1e-6,
  VR7: 50e3, RON_4066: 250, R119: 47e3, C36: 0.033e-6, R120: 22e3,
  VR3: 50e3, R47: 10e3, TM3: 500e3 * 0.318084716796875, R73: 100e3,
  // MEG
  C62: 1e-6, R138: 68e3, VR6: 1e6, R150_R152: 22 + 100, R117: 22e3, VREF: 5.333,
  // VCA output: IC15's output node sees (R127 + R129) || R128 to the bias point
  R_OTA: (47e3 + 22e3) * 220e3 / (47e3 + 22e3 + 220e3), C38: 1e-6, VR8: 50e3,
  // Tune: VR2 (0 to 5.333 V) through R118 against R106 + TM5 at the exponential converter
  TUNE_OCT_PER_V: (27e3 + 5e3 * 0.05665595364634864) / 100e3,
};

// Device numbers taken from the study's models and DC table, not fitted here.
const DEV = {
  HFE21: 400,                  // Q21 baseline hFE (README section 3)
  BETA18: 240, BETA19: 240,    // 2SC536 in the deck
  BETA11: 240,                 // Q11 (2SC536): its base current comes out of Q10's share of I(R73)
  BETA9: 300, IS9: 1.2e-14, VAF9: 60,  // Q9 (2SA733P in the deck), which holds N$561
  VEB10_REF_I: 8.79e-6,        // calibration.json r73_ua at Cutoff 50 %
  V69: 0.598,                  // Q3 Vbe at about 150 uA, the VCA mirror input that D27 feeds
};

// The fitted constants. fit.mjs refits IT, VEB10, KEXP, MIRROR and G_MIX; the README lists what each was
// fitted to. "reference" means results/reference_runs.json, made by reference_runs.py.
const FIT = {
  IT: 0.4668e-3,      // Q21 tail current: loop margins at full Resonance (calibration.json)
  VEB10: 0.5246,      // Q10 Veb at 8.79 uA: R73 current against the Cutoff knob (calibration.json)
  KEXP: 1.0118,       // Q11/Q10 current ratio: ladder current against R73 current (calibration.json)
  V53_DROP: 1.200,    // MEG Darlington drop: N$53 against C62 (envelopes, MEG max 10.13 to 10.28 V)
  V53_SOFT: 0.060,
  V53_MIN: 0.017,     // N$53 at rest, held up by VR5/R61 against Q39 (envelopes, MEG min)
  V_STRIKE: 11.85,    // strike source: C62 peaks 11.33 / 11.43 / 11.48 V at Decay accent / 50 % / 100 %
  T_STRIKE: 1.3e-3,   // strike length: gate to peak 1.2 to 1.4 ms (envelopes_summary.json)
  // VCA control current (reference, Accent 0): a step 5.63 ms after the first gate, 4.74 ms
  // 187.5 ms after the last gate fell and 2.47 ms 62.5 ms after it; 94 % at once, the rest
  // with 3.1 ms; 62.1 uA at 20 ms, decaying with 2.29 s; after the gate, half at 7.7 ms.
  // The step itself grows exponentially, doubling about every 0.13 ms, for 0.7 ms (a whole-
  // chain run sampled every 10 us), and half the step falls at the delays above.
  VEG_I0: 62.65e-6, VEG_TAU: 2.29,
  VEG_DELAY: 5.63e-3, VEG_DELAY_RECOVER: 5.96e-3, VEG_RECOVER_TAU: 0.1, VEG_DELAY_MIN: 0.3e-3,
  VEG_STEP: 0.942, VEG_RISE_TAU: 3.1e-3, VEG_RAMP: 0.70e-3, VEG_RAMP_TAU: 0.19e-3,
  VEG_OFF_HALF: 7.75e-3, VEG_OFF_POW: 5,
  MIRROR: 1.221,      // Q1-Q4 mirror gain with Early effect: accent peaks 210.8 and 115.2 uA (reference)
  G_MIX: 1.557,       // mixer gain to the line out: line Vpp of the reference runs
  // Square from Q8 (reference, four pitches): saw level at the rising edge against pitch,
  // low level, top that droops toward SQ_VINF, and a rising edge 1.4 % of a period long.
  SQ_LOW: 5.375, SQ_START: 8.90, SQ_VINF: 7.07, SQ_TAU: 8.83e-3, SQ_RISE_PHASE: 0.014,
  SQ_EDGE: [[65.41, 8.492], [130.81, 8.251], [261.63, 8.015], [523.25, 7.712]],
};

const SAW_HI = 11.25, SAW_LO = 5.75;  // ideal saw of the study: 5.5 Vpp on 8.5 V, falling ramp
const I117 = (PART.VREF - 0.512) / PART.R117;  // bias of D26, the bottom of the MEG's diode stack
// Line-output volts that make digital full scale. The study's demo pattern peaks near -10 dBFS;
// Volume 10 with Accent 10 goes past it and into the soft limit.
const LINE_FS = 6;

/** Audio taper assumed by the study: t (20 %) of the track at half rotation. */
function taperA(x, t = 0.2) {
  x = Math.min(Math.max(x, 0), 1);
  return x <= 0.5 ? 2 * t * x : t + 2 * (1 - t) * (x - 0.5);
}

/** The two halves of a pot, as the deck makes them (never below 1 ohm). */
function potHalves(total, frac) {
  frac = Math.min(Math.max(frac, 0), 1);
  return [Math.max(total * frac, 1), Math.max(total * (1 - frac), 1)];
}

/**
 * Current through a diode in series with Rs, for a voltage U across both. Returns the
 * current; leaves its slope dI/dU in DS.g and the junction voltage over nvt in DS.w, which
 * a caller can pass back as w0 to start the next solve near the answer.
 */
const DS = { g: 0, w: 0 };
function diodeSeries(U, Rs, nvt, w0) {
  if (U <= 0) {
    const e = Math.exp(Math.max(U / nvt, -60));
    DS.g = IS_D * e / nvt;
    DS.w = Math.log(e);
    return IS_D * (e - 1);
  }
  const wMax = Math.min(U / nvt, Math.log1p(U / (Rs * IS_D)));  // both lie right of the root
  let w = w0 !== undefined && w0 > 0 && w0 < wMax ? w0 : wMax;
  for (let k = 0; k < 50; k++) {
    const ew = Math.exp(w);
    const dw = (Rs * IS_D * (ew - 1) + nvt * w - U) / (Rs * IS_D * ew + nvt);
    w = Math.min(w - dw, wMax);
    if (Math.abs(dw) < 1e-10) break;
  }
  const ew = Math.exp(w);
  DS.g = 1 / (Rs + nvt / (IS_D * ew));
  DS.w = w;
  return IS_D * (ew - 1);
}

function polyBlep(t, dt) {
  if (t < dt) { const x = t / dt; return x + x - x * x - 1; }
  if (t > 1 - dt) { const x = (t - 1) / dt; return x * x + x + x + 1; }
  return 0;
}

// ---------- small dense linear algebra, real and complex ----------
function solveReal(A, b) {  // A: array of rows (copied), b: array; partial pivoting
  const n = b.length;
  const M = A.map((r, i) => r.concat([b[i]]));
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    [M[c], M[p]] = [M[p], M[c]];
    for (let r = c + 1; r < n; r++) {
      const m = M[r][c] / M[c][c];
      for (let k = c; k <= n; k++) M[r][k] -= m * M[c][k];
    }
  }
  const x = new Array(n).fill(0);
  for (let r = n - 1; r >= 0; r--) {
    let s = M[r][n];
    for (let k = r + 1; k < n; k++) s -= M[r][k] * x[k];
    x[r] = s / M[r][r];
  }
  return x;
}

function invReal(A) {
  const n = A.length;
  const cols = [];
  for (let j = 0; j < n; j++) cols.push(solveReal(A, A.map((_, i) => (i === j ? 1 : 0))));
  return A.map((_, i) => cols.map(c => c[i]));
}

const cx = {
  add: (a, b) => [a[0] + b[0], a[1] + b[1]],
  sub: (a, b) => [a[0] - b[0], a[1] - b[1]],
  mul: (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]],
  div: (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; },
  scale: (a, k) => [a[0] * k, a[1] * k],
  abs: a => Math.hypot(a[0], a[1]),
};

function solveComplex(A, b) {
  const n = b.length;
  const M = A.map((r, i) => r.map(z => z.slice()).concat([b[i].slice()]));
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (cx.abs(M[r][c]) > cx.abs(M[p][c])) p = r;
    [M[c], M[p]] = [M[p], M[c]];
    for (let r = c + 1; r < n; r++) {
      const m = cx.div(M[r][c], M[c][c]);
      for (let k = c; k <= n; k++) M[r][k] = cx.sub(M[r][k], cx.mul(m, M[c][k]));
    }
  }
  const x = new Array(n);
  for (let r = n - 1; r >= 0; r--) {
    let s = M[r][n];
    for (let k = r + 1; k < n; k++) s = cx.sub(s, cx.mul(M[r][k], x[k]));
    x[r] = cx.div(s, M[r][r]);
  }
  return x;
}

// ---------- the feedback and output-tap network around the Resonance pot ----------
/**
 * C14, C15, C21, C20 with VR4a, R67, Q18's input, R121, R122 and R124, driven by Q19's
 * emitter (N$144), and Q18's base bias: R99 to N$132, which is C22 || R66 and so no ground
 * below a few hertz. States: C14, C15, C21, C20, C22 voltages. Returns continuous A (5*5),
 * B (5*1), C (2*5), D (2*1) with outputs [N$131 (Q18 base), N$623 (VCA input)].
 */
function tapNetwork(res, P = PART) {
  const [Rlo, Rup] = potHalves(P.VR4, res);  // wiper to ground, wiper to N$27
  const Rq18 = DEV.BETA18 * (P.R94 * (P.R97 + P.R98) / (P.R94 + P.R97 + P.R98));  // Q18's input
  const sq = 1 / P.R67 + 1 / P.R99 + 1 / Rq18;
  const al = (1 / P.R67) / sq, be = (1 / P.R99) / sq;   // n131 = al * n133 + be * n132
  const sig = 1 / P.R121 + 1 / P.R122 + 1 / P.R124;
  const b1 = (1 / P.R121) / sig, b2 = (1 / P.R122) / sig;
  function evalNet(x, u) {
    const n27 = u - x[0];
    const n26 = (n27 / Rup + (x[1] * (1 - al) + be * x[4]) / P.R67 + (x[3] * (1 - b2) + b1 * (n27 - x[2])) / P.R122) /
      (1 / Rup + 1 / Rlo + (1 - al) / P.R67 + (1 - b2) / P.R122);
    const n133 = n26 - x[1];
    const n131 = al * n133 + be * x[4];
    const n25 = n27 - x[2], n21 = n26 - x[3];
    const n623 = b1 * n25 + b2 * n21;
    const iC21 = (n25 - n623) / P.R121, iC20 = (n21 - n623) / P.R122;
    const iC15 = (n133 - n131) / P.R67, iC14 = (n27 - n26) / Rup + iC21;
    const iC22 = (n131 - x[4]) / P.R99 - x[4] / P.R66;
    return [[iC14 / P.C14, iC15 / P.C15, iC21 / P.C21, iC20 / P.C20, iC22 / P.C22], [n131, n623]];
  }
  const n = 5;
  const A = [], C = [[], []];
  for (let i = 0; i < n; i++) A.push([]);
  for (let j = 0; j < n; j++) {
    const x = new Array(n).fill(0);
    x[j] = 1;
    const [dx, y] = evalNet(x, 0);
    for (let i = 0; i < n; i++) A[i][j] = dx[i];
    C[0][j] = y[0];
    C[1][j] = y[1];
  }
  const [dB, yD] = evalNet(new Array(n).fill(0), 1);
  return { A, B: dB.map(v => [v]), C, D: yD.map(v => [v]) };
}

/**
 * Q12's two bases: the VCO through R62 and C17 to N$124, the feedback from Q18's emitter
 * through C23 and R97 to N$60, both returning through R70 and R98 to the bias node N$125,
 * which is C16 || R68 and so no ground below about 7 Hz. States: C17, C16, C23 voltages.
 * Inputs: [N$126, Q18 emitter]. Output: Vd = V(N$124) - V(N$60), Q12's differential input.
 */
function inputNetwork(P, R62, R97) {
  function evalNet(x, u) {
    const n125 = x[1], n128 = u[1] - x[2];
    const n60 = (n128 / R97 + n125 / P.R98) / (1 / R97 + 1 / P.R98);
    const n124 = ((u[0] - x[0]) / R62 + n125 / P.R70) / (1 / R62 + 1 / P.R70);
    const i17 = (n124 - n125) / P.R70, i23 = (n128 - n60) / R97;
    const i16 = i17 + (n60 - n125) / P.R98 - n125 / P.R68;
    return [[i17 / P.C17, i16 / P.C16, i23 / P.C23], [n124 - n60]];
  }
  const A = [[], [], []], B = [[], [], []], C = [[]], D = [[]];
  for (let j = 0; j < 3; j++) {
    const x = [0, 0, 0];
    x[j] = 1;
    const [dx, y] = evalNet(x, [0, 0]);
    for (let i = 0; i < 3; i++) A[i][j] = dx[i];
    C[0][j] = y[0];
  }
  for (let j = 0; j < 2; j++) {
    const u = [0, 0];
    u[j] = 1;
    const [dx, y] = evalNet([0, 0, 0], u);
    for (let i = 0; i < 3; i++) B[i][j] = dx[i];
    D[0][j] = y[0];
  }
  return { A, B, C, D };
}

/** Bilinear discretisation of any state-space (A n*n, B n*m, C p*n, D p*m), flattened rows. */
function discretiseN(net, h) {
  const n = net.A.length, m = net.B[0].length;
  const L = net.A.map((row, i) => row.map((a, j) => (i === j ? 1 : 0) - 0.5 * h * a));
  const R = net.A.map((row, i) => row.map((a, j) => (i === j ? 1 : 0) + 0.5 * h * a));
  const Mi = invReal(L);
  const mul = (X, Y) => X.map(row => Y[0].map((_, j) => row.reduce((acc, v, k) => acc + v * Y[k][j], 0)));
  const Ad = mul(Mi, R);
  const Bd = mul(Mi, net.B.map(row => row.map(b => 0.5 * h * b)));
  const K = mul(net.C, Bd).map((row, i) => row.map((v, j) => v + net.D[i][j]));
  return { n, m, Ad: Float64Array.from(Ad.flat()), Bd: Float64Array.from(Bd.flat()), C: Float64Array.from(net.C.flat()), K: Float64Array.from(K.flat()) };
}

/** Continuous transfer C (sI - A)^-1 B + D at s = j 2 pi f, as complex [p][m]. */
function transfer(net, f) {
  const SA = net.A.map((row, i) => row.map((a, j) => (i === j ? [-a, 2 * Math.PI * f] : [-a, 0])));
  return net.C.map((crow, o) => net.B[0].map((_, j) => {
    const xs = solveComplex(SA, net.B.map(b => [b[j], 0]));
    return xs.reduce((acc, xi, i) => cx.add(acc, cx.scale(xi, crow[i])), [net.D[o][j], 0]);
  }));
}

// ---------- the voice ----------
class X0xVoice {
  /**
   * @param {number} fs host sample rate
   * @param {object} opt os (oversampling), d27 (VCA accent path loads the accent network),
   *   fbMult (feedback multiplier as the study's R97 change), hfe21, inputScale (R62 change),
   *   megHold (hold C62 at a voltage, as the study's AC runs do), d24 (filter accent path),
   *   part (part values replacing PART), fit (constants replacing FIT)
   */
  constructor(fs, opt = {}) {
    this.fs = fs;
    this.os = opt.os || (fs >= 88200 ? 2 : 4);
    this.h = 1 / (fs * this.os);
    this.hHost = 1 / fs;
    this.opt = { d27: true, d24: true, fbMult: 1, hfe21: DEV.HFE21, inputScale: 1, megHold: null, cutTaper: 0.2, ...opt };
    this.knob = { cut: 0.5, res: 0.5, env: 0.5, decay: 0.5, acc: 0.5, tune: 0.5, vol: 0.75, wave: 'saw' };
    this.fit = { ...FIT, ...(opt.fit || {}) };
    this.P = { ...PART, ...(opt.part || {}) };
    this._setupConstants();
    this.reset();
  }

  /**
   * Change the parts the page offers to try, while playing: Q21's hFE, R97 (the feedback
   * divider with R98), and the Cutoff pot's taper (its track fraction at half rotation).
   */
  setCircuit({ hfe21, r97, cutTaper } = {}) {
    if (hfe21) this.opt.hfe21 = hfe21;
    if (r97) this.P.R97 = r97;
    if (cutTaper && cutTaper !== this.opt.cutTaper) {
      // A builder trims TM3 for the kit's 500 Hz at Cutoff 50 %, so with another taper TM3 is
      // retrimmed to give the same R73 current there; the top of the knob then moves.
      const P = this.P, V147 = PART.V561 + this.fit.VEB10;
      const at = (t) => {
        const [Ras, Rse] = potHalves(P.VR3, t), r3 = Ras + P.R47, Rp = 1 / (1 / r3 + 1 / Rse);
        return [(PART.V561 / r3 + 12 / Rse) * Rp, Rp];
      };
      const [V0, Rp0] = at(this.opt.cutTaper), [V1, Rp1] = at(cutTaper);
      P.TM3 = (V1 - V147) * (Rp0 + P.TM3 + P.R73) / (V0 - V147) - Rp1 - P.R73;
      this.opt.cutTaper = cutTaper;
    }
    this._setupLoop();
    this._updateKnobs();
  }

  _setupConstants() {
    this._setupLoop();
    this._setupOutput();
  }

  /** The feedback loop's constants. The states stay: they are capacitor voltages. */
  _setupLoop() {
    const P = this.P;
    const o = this.opt;
    this.R97 = (P.R97 + P.R98) / o.fbMult - P.R98;
    this.Rin = 1 / (1 / P.R114 + this.fit.IT / (2 * o.hfe21 * VT));  // per side, Q21 base || R114
    this.A21h = P.R110 * this.fit.IT / 2;                            // Q21A collector swing
    this.g19 = P.R96 / (P.R96 + VT / (5.56 / P.R96));          // emitter followers
    this.g18 = 1 / (1 + VT / (4.36 / P.R94) / (P.R94 * (this.R97 + P.R98) / (P.R94 + this.R97 + P.R98)));
    this.R62 = P.R62 / o.inputScale;
    const g = (tau) => { const gg = this.h / (2 * tau); return gg / (1 + gg); };
    this.G25 = g(this.Rin * P.C25);
    this.inNet = inputNetwork(P, this.R62, this.R97);
    this.inD = discretiseN(this.inNet, this.h);
    this._netRes = -1;
  }

  _setupOutput() {
    const P = this.P;
    const gh = this.hHost / (2 * P.C38 * P.VR8);
    this.G38 = gh / (1 + gh);
    // decimation FIR: windowed sinc at the oversampled rate, cutoff 0.45 * fs host
    const N = this.os === 1 ? 1 : 16 * this.os + 1;
    this.fir = new Float64Array(N);
    if (N === 1) this.fir[0] = 1;
    else {
      let s = 0;
      const fc = 0.45 / this.os;
      for (let i = 0; i < N; i++) {
        const m = i - (N - 1) / 2;
        const sinc = m === 0 ? 2 * fc : Math.sin(2 * Math.PI * fc * m) / (Math.PI * m);
        const w = 0.42 - 0.5 * Math.cos(2 * Math.PI * i / (N - 1)) + 0.08 * Math.cos(4 * Math.PI * i / (N - 1));
        this.fir[i] = sinc * w;
        s += this.fir[i];
      }
      for (let i = 0; i < N; i++) this.fir[i] /= s;
    }
    this.ring = new Float64Array(N);
    this.ringPos = 0;
  }

  reset() {
    const P = this.P;
    // ladder and loop
    this.v = new Float64Array(4);
    this.fo = new Float64Array(4);
    this.xi = new Float64Array(3); this.zi = new Float64Array(3);   // C17, C16, C23
    this.primeInput(8.5);      // C17 charged to the saw's mean, as with the VCO running
    this.s25 = 0;
    this.xn = new Float64Array(5); this.zn = new Float64Array(5); this.u144 = 0;
    this.out623 = 0; this.out144 = 0; this.vd = 0;
    // control, absolute volts
    this.V46 = 0.178;          // C62 at the deck's operating point, before any note
    this.iDis = 0; this.wDis = undefined; this.wCh = undefined; this.w24 = undefined; this.w27 = undefined;
    this.V161 = P.V561;     // C13 settles at the summing node
    this.V34 = 0;              // C36
    this.V35 = 0; this.V152 = P.V561; this.i24 = 0; this.i27 = 0;
    this.vsum = 0; this.V561 = P.V561; this.V134 = P.V561; this.IR73 = 8e-6;
    this.IL = 7e-6; this.ILprev = 7e-6;
    // lines and timers
    this.gate = 0; this.accent = 0; this.slide = 0;
    this.strikeLeft = 0;
    this.time = 0; this.lastOff = -1e9;
    this.veg = { tOn: -1e9, tOff: -1e9, start: 1e9 };
    this.vegPrev = { tOn: -1e9, tOff: -1e9, start: 1e9 };
    this.iabc = 0;
    this.cvDac = 2.0; this.cv = 2.0; this.freqOverride = 0;
    this.phase = 0; this.sqT = 0; this.sqHigh = false; this.sqLevel = FIT.SQ_LOW; this.sqPrev = FIT.SQ_LOW;
    this.s38 = 0;
    this.events = [];
    this.sampleIndex = 0;
    this.extInput = null;
    this.ring.fill(0);
    this._controlInit();
  }

  /** Settle the control network with no notes, as the deck's operating point does. */
  _controlInit() {
    this._updateKnobs();
    for (let i = 0; i < 2000; i++) this._controlStep();
    // C13 (about 0.15 s) and C36 sit at their DC points in the deck's operating point
    this.V161 = this.V152 = this.V134;
    for (let i = 0; i < 2000; i++) this._controlStep();
    this.ILprev = this.IL;
  }

  // ---------- events (firmware lines) ----------
  /** Schedule a line change at time t seconds from the start of the render. */
  schedule(t, type, value) {
    const n = Math.round(t * this.fs);
    let i = this.events.length;
    while (i > 0 && this.events[i - 1].n > n) i--;
    this.events.splice(i, 0, { n, type, value });
  }

  _apply(type, value) {
    if (type === 'gate') {
      if (value && !this.gate) {         // a struck note; a tied one keeps the gate high
        const f = this.fit, t = this.time;
        this.strikeLeft = f.T_STRIKE;
        const delay = Math.max(f.VEG_DELAY_MIN, f.VEG_DELAY - f.VEG_DELAY_RECOVER * Math.exp(-(t - this.lastOff) / f.VEG_RECOVER_TAU));
        this.vegPrev = this.veg;
        // the ramp starts so that it is half way at the delay
        this.veg = { tOn: t, tOff: 1e9, start: t + delay - (f.VEG_RAMP + f.VEG_RAMP_TAU * Math.log(0.5)) };
      }
      if (!value && this.gate) { this.veg.tOff = this.time; this.lastOff = this.time; }
      this.gate = value ? 1 : 0;
    } else if (type === 'accent') this.accent = value ? 1 : 0;
    else if (type === 'slide') this.slide = value ? 1 : 0;
    else if (type === 'cv') { this.cvDac = value; if (!this.slide) this.cv = value; }
    else if (type === 'freq') this.freqOverride = value;
    else if (type in this.knob) { this.knob[type] = value; this._updateKnobs(); }
  }

  // ---------- control path, host rate ----------
  /** Everything in the control path that depends only on the knobs. */
  _updateKnobs() {
    const P = this.P;
    const k = this.knob, h = this.hHost, c = this.kc || (this.kc = {});
    c.RdisNorm = P.R138 + P.VR6 * taperA(k.decay);
    c.RdisAcc = P.R138 + P.RON_4066;
    const [Ra, Rb] = potHalves(P.VR4, k.res);         // VR4b: N$160 to wiper, wiper to C13
    const [Rlo7, Rup7] = potHalves(P.VR7, k.acc);
    const ru = Rup7 + P.RON_4066;
    c.Rth7on = Rlo7 * ru / (Rlo7 + ru);
    c.k7on = Rlo7 / (Rlo7 + ru);
    c.Rth7off = Rlo7;
    c.Rs1 = P.R46 + Ra;
    c.Rbd = Rb + h / P.C13;
    c.R36 = h / P.C36;
    const [Ras5, Rse5] = potHalves(P.VR5, taperA(k.env));
    c.r1 = Ras5 + P.R61;
    c.Rse5 = Rse5;
    c.Rp5 = 1 / (1 / c.r1 + 1 / Rse5);
    c.Rt5 = c.Rp5 + P.R63;
    const [Ras3, Rse3] = potHalves(P.VR3, taperA(k.cut, this.opt.cutTaper));
    c.r3 = Ras3 + P.R47;
    c.Rse3 = Rse3;
    c.Rp3 = 1 / (1 / c.r3 + 1 / Rse3);
    c.Rt3 = c.Rp3 + P.TM3 + P.R73;
  }

  _controlStep() {
    const P = this.P;
    const h = this.hHost, f = this.fit, kc = this.kc;
    // MEG: C62 struck through D37 and R150 + R152, discharged through R138 (+ VR6) into D28, D29, D26
    const Rdis = this.accent ? kc.RdisAcc : kc.RdisNorm;
    if (this.opt.megHold !== null) {
      this.V46 = this.opt.megHold;
    } else {
      const v50 = NVT_D * Math.log1p((I117 + this.iDis) / IS_D);
      this.iDis = Math.max(diodeSeries(this.V46 - v50, Rdis, 2 * NVT_D, this.wDis), 0);
      this.wDis = DS.w;
      let iCh = 0;
      if (this.strikeLeft > 0) {
        const kh = h / P.C62;
        iCh = Math.max(diodeSeries(f.V_STRIKE - this.V46 + kh * this.iDis, P.R150_R152 + kh, NVT_D, this.wCh), 0);
        this.wCh = DS.w;
        this.strikeLeft -= h;
      }
      this.V46 += h / P.C62 * (iCh - this.iDis);
    }
    const x = (this.V46 - f.V53_DROP) / f.V53_SOFT;
    const V53 = f.V53_MIN + f.V53_SOFT * (x > 30 ? x : Math.log1p(Math.exp(x)));
    this.V53 = V53;
    // Accent network, backward Euler, unknowns V35 (a), V152 (c), V34 (e)
    const Rth7 = this.accent ? kc.Rth7on : kc.Rth7off;
    const Vth7 = this.accent ? V53 * kc.k7on : 0;
    const Rs1 = kc.Rs1, Rbd = kc.Rbd, R36 = kc.R36;
    const d24 = this.opt.d24, d27 = this.opt.d27;
    let a = this.V35, c = this.V152, e = this.V34;
    const V134 = this.V134;
    let i24 = 0, g24 = 0, i27 = 0, g27 = 0;
    for (let it = 0; it < 30; it++) {
      if (d24) { i24 = diodeSeries(a - c, Rs1, NVT_D, this.w24); g24 = DS.g; this.w24 = DS.w; }
      if (d27) { i27 = diodeSeries(e - DEV.V69, P.R120, NVT_D, this.w27); g27 = DS.g; this.w27 = DS.w; }
      const Fa = (Vth7 - a) / Rth7 + (e - a) / P.R119 - i24;
      const Fc = i24 + (V134 - c) / P.R72 + (this.V161 - c) / Rbd;
      const Fe = (a - e) / P.R119 + (this.V34 - e) / R36 - i27;
      if (Math.abs(Fa) + Math.abs(Fc) + Math.abs(Fe) < 1e-12) break;
      const Jaa = -1 / Rth7 - 1 / P.R119 - g24, Jac = g24, Jae = 1 / P.R119;
      const Jca = g24, Jcc = -g24 - 1 / P.R72 - 1 / Rbd;
      const Jea = 1 / P.R119, Jee = -1 / P.R119 - 1 / R36 - g27;
      // eliminate c and e into the a row
      const da = (-Fa + Jac * Fc / Jcc + Jae * Fe / Jee) / (Jaa - Jac * Jca / Jcc - Jae * Jea / Jee);
      a += da; c += (-Fc - Jca * da) / Jcc; e += (-Fe - Jea * da) / Jee;
    }
    this.V161 += (c - this.V161) / Rbd * (h / P.C13);
    this.V35 = a; this.V152 = c; this.V34 = e; this.i24 = i24; this.i27 = Math.max(i27, 0);
    const iR72 = (c - V134) / P.R72;
    // N$561 is Q9's emitter. Everything below is referred to it, and it moves with Q9's current.
    const V561 = this.V561;
    // Env Mod: VR5 between N$49 (R61 to N$561) and the MEG, wiper through R63
    const Vth5 = (V561 / kc.r1 + V53 / kc.Rse5) * kc.Rp5;
    // summing node: Vs = R71 * (I(R63) + I(R72)), with I(R63) depending on Vs itself
    const vs = P.R71 * ((Vth5 - V561) / kc.Rt5 + iR72) / (1 + P.R71 / kc.Rt5);
    this.vsum = vs;
    this.V134 = V561 + vs;
    // Cutoff: VR3 between N$142 (R47 to N$561) and +12 V, wiper through TM3 and R73 to Q10
    const Vth3 = (V561 / kc.r3 + 12 / kc.Rse3) * kc.Rp3;
    let I = this.IR73;
    for (let it = 0; it < 2; it++) {
      const V147 = this.V134 + f.VEB10 + VT * Math.log(Math.max(I, 1e-9) / DEV.VEB10_REF_I);
      I = (Vth3 - V147) / kc.Rt3;
    }
    this.IR73 = I;
    // Q11's collector current is K * I(Q10) * exp(Vs/VT), and I(Q10) = I(R73) - Q11's base current
    const ke = f.KEXP * Math.exp(vs / VT);
    this.IL = Math.max(ke * I / (1 + ke / DEV.BETA11), 1e-9);
    // Q9 sinks the VR3 divider, the VR5 path, R71, Q10 and Q11: its Veb, and so N$561, follow
    const iR63 = (Vth5 - this.V134) / kc.Rt5;
    const i61 = (Vth5 - iR63 * kc.Rp5 - V561) / kc.r1;
    const i47 = (Vth3 - I * kc.Rp3 - V561) / kc.r3;
    const iE9 = Math.max(i47 + i61 + vs / P.R71 + I + this.IL, 1e-6);
    this.V561 = P.VREF / 2 + P.R64_R65 * iE9 / (DEV.BETA9 + 1) +
      VT * Math.log(iE9 * DEV.BETA9 / (DEV.BETA9 + 1) / (DEV.IS9 * (1 + V561 / DEV.VAF9)));
    // VCA control current: the VEG with the reference timing, plus D27's accent current.
    // A note struck while the last one is still closing keeps that tail until its own step.
    const t = this.time;
    const veg = Math.max(this._vegLevel(this.veg, t), t < this.veg.start ? this._vegLevel(this.vegPrev, t) : 0);
    this.iabc = veg + (this.opt.d27 ? f.MIRROR * this.i27 : 0);
    // slide: the CV reaches the VCO through 22 ms while SLIDE is high
    if (this.slide) this.cv += (this.cvDac - this.cv) * (1 - Math.exp(-h / 0.022));
    else this.cv = this.cvDac;
    this.time += h;
  }

  _vegLevel(n, t) {
    if (t < n.start) return 0;
    const f = this.fit, ts = t - n.start;
    const tr = ts - f.VEG_RAMP;
    const on = tr < 0 ? f.VEG_STEP * Math.exp(tr / f.VEG_RAMP_TAU) : 1 - (1 - f.VEG_STEP) * Math.exp(-tr / f.VEG_RISE_TAU);
    const off = t > n.tOff ? 1 / (1 + Math.pow((t - n.tOff) / f.VEG_OFF_HALF, f.VEG_OFF_POW)) : 1;
    return f.VEG_I0 * on * off * Math.exp(-(t - n.tOn) / f.VEG_TAU);
  }

  /** Re-discretise the Resonance-pot network when the knob has moved. */
  _updateNet() {
    const r = this.knob.res;
    if (Math.abs(r - this._netRes) < 1e-4) return;
    this._netRes = r;
    this.net = discretiseN(tapNetwork(r, this.P), this.h);
  }

  // ---------- oscillator, oversampled rate ----------
  _oscStep() {
    const P = this.P;
    if (this.extInput) return this.extInput();
    const f = this.freqOverride || 65.41 * Math.pow(2, this.cv - 2 + (this.knob.tune - 0.5) * P.VREF * P.TUNE_OCT_PER_V);
    const dt = f * this.h;
    this.phase += dt;
    if (this.phase >= 1) this.phase -= 1;
    const p = this.phase;
    if (this.knob.wave !== 'square') return 8.5 + 2.75 * (1 - 2 * p) + 2.75 * polyBlep(p, dt);
    // Square from Q8: low until the saw falls past a pitch-dependent level, then a top that
    // droops toward SQ_VINF, falling back when the saw resets.
    const F = this.fit;
    if (f !== this.sqF) {   // saw level at the rising edge, interpolated in octaves
      const E = F.SQ_EDGE, L = Math.log2(f);
      let j = 0;
      while (j < E.length - 2 && L > Math.log2(E[j + 1][0])) j++;
      const l0 = Math.log2(E[j][0]), l1 = Math.log2(E[j + 1][0]);
      const vth = E[j][1] + (E[j + 1][1] - E[j][1]) * (L - l0) / (l1 - l0);
      this.sqF = f;
      this.sqPr = Math.min(Math.max((SAW_HI - vth) / (SAW_HI - SAW_LO), 0.05), 0.95);
    }
    const pr = this.sqPr;
    let y;
    if (p >= pr) {
      if (!this.sqHigh) { this.sqHigh = true; this.sqT = (p - pr) / f; }
      else this.sqT += this.h;
      const top = F.SQ_VINF + (F.SQ_START - F.SQ_VINF) * Math.exp(-this.sqT / F.SQ_TAU);
      const ramp = Math.min(this.sqT * f / F.SQ_RISE_PHASE, 1);
      this.sqLevel = top;
      y = F.SQ_LOW + (top - F.SQ_LOW) * ramp * ramp * (3 - 2 * ramp);
    } else {
      this.sqHigh = false;
      y = F.SQ_LOW;
    }
    // the fall at the saw reset is fast: band-limit it with a BLEP of the drop
    const drop = this.sqLevel - F.SQ_LOW;
    if (p < dt) y += -0.5 * drop * polyBlep(p, dt);
    else if (p > 1 - dt && p >= pr) y += -0.5 * drop * polyBlep(p, dt);
    return y;
  }

  // ---------- the ladder and its loop, oversampled rate ----------
  _coreStep(u126, I) {
    const P = this.P;
    const h = this.h, hh = 0.5 * h;
    // affine maps of the loop for this step
    const a0 = 1 - this.G25, b0 = -(1 - this.G25) * this.s25;
    const N = this.net, x = this.xn, z = this.zn, Ad = N.Ad, Bd = N.Bd, NC = N.C, nn = N.n;
    let bA = 0, bB = 0;
    for (let i = 0; i < nn; i++) {
      let acc = Bd[i] * this.u144;
      for (let j = 0; j < nn; j++) acc += Ad[nn * i + j] * x[j];
      z[i] = acc;
      bA += NC[i] * acc;
      bB += NC[nn + i] * acc;
    }
    const k0 = N.K[0], k1n = N.K[1];
    const g19 = this.g19, g18 = this.g18;
    // Q12's input: Vd = vin + Kv * Vc, from the input network with this step's VCO sample
    const Q = this.inD, xi = this.xi, zi = this.zi, ui1 = this.ui1, ui2 = this.ui2;
    for (let i = 0; i < 3; i++) zi[i] = Q.Ad[3 * i] * xi[0] + Q.Ad[3 * i + 1] * xi[1] + Q.Ad[3 * i + 2] * xi[2] + Q.Bd[2 * i] * ui1 + Q.Bd[2 * i + 1] * ui2;
    const vin = Q.C[0] * zi[0] + Q.C[1] * zi[1] + Q.C[2] * zi[2] + Q.K[0] * u126 + Q.K[1] * g18 * bA;
    const Kfb = -Q.K[1] * g18 * k0 * g19, Bfb = 0;
    const A21h = this.A21h, Rin = this.Rin;
    const k1 = 1 / (2 * P.C18), k = 1 / (2 * P.C);
    const v = this.v, fo = this.fo;
    let v1 = v[0], v2 = v[1], v3 = v[2], v4 = v[3];
    let f1 = 0, f2 = 0, f3 = 0, f4 = 0, Vc = 0, Vd = 0;
    for (let it = 0; it < 12; it++) {
      const vb = a0 * v4 + b0;
      const t21 = Math.tanh(vb * INV2VT);
      Vc = -A21h * t21;
      Vd = vin - (Kfb * Vc + Bfb);
      const t0 = Math.tanh(Vd * INV2VT);
      const t1 = Math.tanh((v2 - v1) * INV2VT);
      const t2 = Math.tanh((v3 - v2) * INV2VT);
      const t3 = Math.tanh((v4 - v3) * INV2VT);
      const t4 = Math.tanh(-v4 * INV2VT);
      f1 = (I * t1 - I * t0) * k1;
      f2 = (I * t2 - I * t1) * k;
      f3 = (I * t3 - I * t2) * k;
      f4 = (I * t4 - I * t3 - vb / Rin) * k;
      const r1 = v1 - v[0] - hh * (fo[0] + f1);
      const r2 = v2 - v[1] - hh * (fo[1] + f2);
      const r3 = v3 - v[2] - hh * (fo[2] + f3);
      const r4 = v4 - v[3] - hh * (fo[3] + f4);
      if (Math.abs(r1) + Math.abs(r2) + Math.abs(r3) + Math.abs(r4) < 1e-13) break;
      const d0 = I * (1 - t0 * t0) * INV2VT, d1 = I * (1 - t1 * t1) * INV2VT, d2 = I * (1 - t2 * t2) * INV2VT;
      const d3 = I * (1 - t3 * t3) * INV2VT, d4 = I * (1 - t4 * t4) * INV2VT;
      const dVd4 = Kfb * A21h * (1 - t21 * t21) * INV2VT * a0;  // dVd/dv4
      // Jacobian of the residual: tridiagonal plus the loop term in row 1, column 4
      const J11 = 1 + hh * k1 * d1, J12 = -hh * k1 * d1, J14 = hh * k1 * d0 * dVd4;
      const J21 = -hh * k * d1, J22 = 1 + hh * k * (d1 + d2), J23 = -hh * k * d2;
      const J32 = -hh * k * d2, J33 = 1 + hh * k * (d2 + d3), J34 = -hh * k * d3;
      const J43 = -hh * k * d3, J44 = 1 + hh * k * (d3 + d4 + a0 / Rin);
      // elimination: row 2 by row 1, row 3 by row 2, row 4 by row 3 (fill-in in column 4)
      const m2 = J21 / J11;
      const B22 = J22 - m2 * J12, B24 = -m2 * J14, q2 = -r2 + m2 * r1;
      const m3 = J32 / B22;
      const B33 = J33 - m3 * J23, B34 = J34 - m3 * B24, q3 = -r3 - m3 * q2;
      const m4 = J43 / B33;
      const B44 = J44 - m4 * B34, q4 = -r4 - m4 * q3;
      const x4 = q4 / B44;
      const x3 = (q3 - B34 * x4) / B33;
      const x2 = (q2 - J23 * x3 - B24 * x4) / B22;
      const x1 = (-r1 - J12 * x2 - J14 * x4) / J11;
      v1 += x1; v2 += x2; v3 += x3; v4 += x4;
    }
    v[0] = v1; v[1] = v2; v[2] = v3; v[3] = v4;
    fo[0] = f1; fo[1] = f2; fo[2] = f3; fo[3] = f4;
    // advance the linear parts with the solved values
    const vb25 = (v4 - this.s25) * this.G25;
    const lp25 = vb25 + this.s25;
    this.s25 = lp25 + vb25;
    const u144 = g19 * Vc;
    for (let i = 0; i < nn; i++) x[i] = z[i] + Bd[i] * u144;
    this.u144 = u144;
    const n131 = k0 * u144 + bA;
    this.out623 = k1n * u144 + bB;
    const v130 = g18 * n131;
    for (let i = 0; i < 3; i++) xi[i] = zi[i] + Q.Bd[2 * i] * u126 + Q.Bd[2 * i + 1] * v130;
    this.ui1 = u126; this.ui2 = v130;
    this.out144 = u144;
    this.vd = Vd;
  }

  /**
   * Render n host samples into out (Float32Array or Float64Array). probes: optional object of
   * arrays filled per host sample: ua (ladder current, uA), vsum (N$134 - N$561, mV), v53
   * and v46 (MEG output and C62, V), v623 (VCA input), v144 (Q19 emitter), vd (Q12 input),
   * src (VCO), iabc (VCA control current, uA), v161 (C13), v34 (C36), line and v48 (V).
   */
  process(out, n, probes) {
    const P = this.P;
    this._updateNet();
    this._updateKnobs();
    const os = this.os, fir = this.fir, ring = this.ring, NF = fir.length;
    const volG = taperA(this.knob.vol) * this.fit.G_MIX / LINE_FS;
    for (let i = 0; i < n; i++) {
      while (this.events.length && this.events[0].n <= this.sampleIndex) {
        const ev = this.events.shift();
        this._apply(ev.type, ev.value);
        if (ev.type === 'res') this._updateNet();
      }
      this.ILprev = this.IL;
      this._controlStep();
      const I0 = this.ILprev, I1 = this.IL, ia = this.iabc;
      for (let s = 1; s <= os; s++) {
        const I = I0 + (I1 - I0) * (s / os);
        const u = this._oscStep();
        this.src = u;
        this._coreStep(u, I);
        // IC15 as the study's OTA: Iabc * tanh(Vin / 2VT) into the output node
        ring[this.ringPos] = ia * P.R_OTA * Math.tanh(this.out623 * INV2VT);
        this.ringPos = (this.ringPos + 1) % NF;
      }
      let acc = 0;
      for (let j = 0, p = this.ringPos; j < NF; j++) {
        acc += fir[j] * ring[p];
        p = p + 1 === NF ? 0 : p + 1;
      }
      // C38 into the volume pot, then the mixer
      const vv = (acc - this.s38) * this.G38;
      const lp = vv + this.s38;
      this.s38 = lp + vv;
      const y = (acc - lp) * volG;
      // soft limit above 0.7 of full scale: continuous, and with the same slope at the knee
      const ay = Math.abs(y);
      out[i] = ay <= 0.7 ? y : Math.sign(y) * (0.7 + 0.3 * Math.tanh((ay - 0.7) / 0.3));
      if (probes) {
        if (probes.ua) probes.ua[i] = this.IL * 1e6;
        if (probes.vsum) probes.vsum[i] = this.vsum * 1e3;
        if (probes.v53) probes.v53[i] = this.V53;
        if (probes.v46) probes.v46[i] = this.V46;
        if (probes.v623) probes.v623[i] = this.out623;
        if (probes.v144) probes.v144[i] = this.out144;
        if (probes.vd) probes.vd[i] = this.vd;
        if (probes.src) probes.src[i] = this.src;
        if (probes.iabc) probes.iabc[i] = this.iabc * 1e6;
        if (probes.v161) probes.v161[i] = this.V161;
        if (probes.v34) probes.v34[i] = this.V34;
        if (probes.line) probes.line[i] = (acc - lp) * taperA(this.knob.vol) * this.fit.G_MIX;
        if (probes.v48) probes.v48[i] = acc;
      }
      this.sampleIndex++;
    }
  }

  /** Ladder current and summing-node voltage now, for the page's meter. */
  get ladderMicroamps() { return this.IL * 1e6; }

  /**
   * Charge C17 to a source level, as the deck's operating point does: the study's runs start
   * with the saw at its top (11.25 V) or the square low; a running VCO sits at its mean.
   */
  primeInput(volts) { this.xi.fill(0); this.xi[0] = volts; this.ui1 = volts; this.ui2 = 0; }

  /**
   * Small-signal response of this model at its present control state, analytically:
   * complex V623/V126, V144/V126 and the loop gain T, broken at Q18's emitter, a low
   * impedance (T = +1 oscillates). Ts is the loop gain broken as the study broke it, between
   * C23 and R97 with R97 + R98 as the return's load, for comparing margins. freqs in Hz.
   */
  smallSignal(freqs) {
    const P = this.P;
    const I = this.IL, g = I / (2 * VT);
    const net = tapNetwork(this.knob.res, this.P);
    const tau25 = this.Rin * P.C25;
    // the study's break: a source drives R97 (C23 shorted), the return is C23 into R97 + R98
    const injNet = inputNetwork({ ...P, C23: 1e3 }, this.R62, this.R97);
    const tauRet = P.C23 * (this.R97 + P.R98);
    const A21 = this.A21h * INV2VT;
    const out = [];
    for (const f of freqs) {
      const s = [0, 2 * Math.PI * f];
      const hp = tau => { const st = cx.scale(s, tau); return cx.div(st, cx.add([1, 0], st)); };
      // ladder with Q21 loading, driven by Vd = 1
      const H25 = hp(tau25);
      const c1 = cx.scale(s, 2 * P.C18), c = cx.scale(s, 2 * P.C);
      const z = [0, 0];
      const M = [
        [cx.add(c1, [g, 0]), [-g, 0], z, z],
        [[-g, 0], cx.add(c, [2 * g, 0]), [-g, 0], z],
        [z, [-g, 0], cx.add(c, [2 * g, 0]), [-g, 0]],
        [z, z, [-g, 0], cx.add(cx.add(c, [2 * g, 0]), cx.scale(H25, 1 / this.Rin))],
      ];
      const lad = solveComplex(M, [[-g, 0], z, z, z]);
      const L = lad[3];
      const Hn = transfer(net, f).map(row => row[0]);           // [N$131, N$623] over N$144
      const fwd = cx.scale(H25, -A21 * this.g19);              // v4 -> V144
      const Hi = transfer(this.inNet, f)[0];                     // [Vd/V126, Vd/V(Q18 emitter)]
      const Fb = cx.scale(cx.mul(cx.mul(fwd, Hn[0]), Hi[1]), this.g18);   // v4 -> Vd through the loop
      const T = cx.mul(L, Fb);
      const Hj = transfer(injNet, f)[0][1];
      const Ts = cx.mul(cx.mul(cx.mul(cx.mul(Hj, L), fwd), Hn[0]), cx.scale(hp(tauRet), this.g18));
      const vin = Hi[0];
      const v4 = cx.div(cx.mul(L, vin), cx.sub([1, 0], T));
      out.push({ f, vca: cx.mul(cx.mul(v4, fwd), Hn[1]), top: cx.mul(v4, fwd), T, Ts });
    }
    return out;
  }
}

// ---------- the sequencer: the stock firmware's note_on / note_off timing ----------
/**
 * Pattern steps: { note: 0..12 semitones above C, oct: -1|0|1, accent, slide, rest }.
 * note_on (step start): SLIDE takes the previous step's slide flag, ACCENT this step's;
 * a note latches the DAC and raises the gate. note_off (half a step later): SLIDE rises if
 * this step slides, otherwise the gate falls. Rests neither latch nor strike.
 */
class Sequencer {
  constructor(voice) {
    this.voice = voice;
    this.steps = [];
    this.bpm = 120;
    this.running = false;
    this.pos = 0;
    this.nextStepSample = 0;
    this.half = false;
    this.prevSlide = false;
    this.current = -1;
  }

  start(atSample) { this.running = true; this.pos = 0; this.half = false; this.nextStepSample = atSample; this.prevSlide = false; }

  stop(atSample) {
    this.running = false;
    const v = this.voice;
    v.events = v.events.filter(e => e.n < atSample);
    v.events.push({ n: atSample, type: 'gate', value: 0 }, { n: atSample, type: 'slide', value: 0 }, { n: atSample, type: 'accent', value: 0 });
    this.current = -1;
  }

  /** One note, as note_on and note_off would play it, while the pattern is stopped. */
  audition(atSample, st) {
    if (this.running || st.rest) return;
    const v = this.voice;
    const half = Math.round(v.fs * 60 / this.bpm / 8);
    v.events = v.events.filter(e => e.n < atSample);
    v.events.push({ n: atSample, type: 'slide', value: 0 }, { n: atSample, type: 'accent', value: st.accent ? 1 : 0 },
      { n: atSample, type: 'cv', value: 2.0 + (st.note + 12 * (st.oct || 0)) / 12 }, { n: atSample + 1, type: 'gate', value: 1 },
      { n: atSample + half, type: 'gate', value: 0 });
    v.events.sort((a, b) => a.n - b.n);
  }

  /** Queue the line changes that fall inside the next n host samples. */
  plan(fromSample, n) {
    if (!this.running || !this.steps.length) return;
    const v = this.voice;
    const stepLen = v.fs * 60 / this.bpm / 4;
    while (this.nextStepSample < fromSample + n) {
      const t = Math.max(Math.round(this.nextStepSample), fromSample);
      const st = this.steps[this.pos % this.steps.length];
      if (!this.half) {
        v.events.push({ n: t, type: 'slide', value: this.prevSlide ? 1 : 0 });
        v.events.push({ n: t, type: 'accent', value: st.accent ? 1 : 0 });
        if (!st.rest) {
          v.events.push({ n: t, type: 'cv', value: 2.0 + (st.note + 12 * (st.oct || 0)) / 12 });
          v.events.push({ n: t + 1, type: 'gate', value: 1 });
        }
        this.current = this.pos % this.steps.length;
        this.half = true;
        this.nextStepSample += stepLen / 2;
      } else {
        if (!st.rest && st.slide) v.events.push({ n: t, type: 'slide', value: 1 });
        else v.events.push({ n: t, type: 'gate', value: 0 });
        this.prevSlide = !st.rest && !!st.slide;
        this.half = false;
        this.pos++;
        this.nextStepSample += stepLen / 2;
      }
    }
    v.events.sort((a, b) => a.n - b.n);
  }
}

const api = { X0xVoice, Sequencer, taperA, PART, DEV, FIT, VT, tapNetwork, inputNetwork, discretiseN, transfer, diodeSeries };

/**
 * Run every path of the code for a moment before the first real block, so the JIT has
 * compiled it: otherwise the first notes, accents and slides can each take longer than a
 * block and the output crackles. Takes about 60 ms on a laptop.
 */
function warmUp(fs, os) {
  const v = new X0xVoice(fs, { os });
  const s = new Sequencer(v);
  s.steps = [{ note: 0, oct: 0, accent: true, slide: false, rest: false }, { note: 12, oct: 0, accent: false, slide: true, rest: false },
    { note: 3, oct: 1, accent: true, slide: true, rest: false }, { note: 0, oct: -1, accent: false, slide: false, rest: true }];
  s.bpm = 240;
  s.start(0);
  const buf = new Float32Array(128);
  const n = Math.ceil(0.6 * fs / 128);
  for (let b = 0; b < n; b++) {
    if (b === n >> 1) v.knob.wave = 'square';
    v.knob.res = 0.5 + 0.5 * Math.sin(b / 9);
    v.knob.cut = 0.5 + 0.5 * Math.sin(b / 7);
    v.knob.env = 0.5 + 0.5 * Math.sin(b / 5);
    s.plan(b * 128, 128);
    v.process(buf, 128);
  }
  s.stop(n * 128);
  s.audition(n * 128 + 64, s.steps[0]);
  for (let b = n; b < n + 40; b++) v.process(buf, 128);
}

const KNOB_NAMES = ['cut', 'res', 'env', 'decay', 'acc', 'tune', 'vol'];
const clock = () => (globalThis.performance && performance.now ? performance.now() : Date.now());

if (typeof registerProcessor === 'function') {
  class X0xProcessor extends AudioWorkletProcessor {
    constructor(options) {
      super();
      const o = (options && options.processorOptions) || {};
      warmUp(sampleRate, o.os);
      this.voice = new X0xVoice(sampleRate, { os: o.os });
      this.budgetMs = 1000 * 128 / sampleRate;
      this.worstMs = 0;
      this.late = 0;
      this.seq = new Sequencer(this.voice);
      this.target = { ...this.voice.knob };
      this.frame = 0;
      this.lastStep = -2;
      this.meterEvery = Math.round(sampleRate / 30);
      this.meterCount = 0;
      this.peakUa = 0;
      this.busyMs = 0;
      // the page's state at creation, so the first block already has it
      if (o.knobs) Object.assign(this.target, o.knobs), Object.assign(this.voice.knob, o.knobs);
      if (o.wave) this.voice.knob.wave = o.wave;
      if (o.steps) this.seq.steps = o.steps;
      if (o.bpm) this.seq.bpm = o.bpm;
      if (o.circuit) this.voice.setCircuit(o.circuit);
      if (o.play) this.seq.start(64);
      this.port.onmessage = (e) => this.onMessage(e.data);
    }

    onMessage(m) {
      const v = this.voice;
      if (m.type === 'knob') {
        if (m.name === 'wave') v.knob.wave = m.value;
        else this.target[m.name] = m.value;
      } else if (m.type === 'pattern') this.seq.steps = m.steps;
      else if (m.type === 'tempo') this.seq.bpm = m.bpm;
      else if (m.type === 'play') this.seq.start(this.frame + 64);
      else if (m.type === 'stop') this.seq.stop(this.frame);
      else if (m.type === 'audition') this.seq.audition(this.frame + 64, m.step);
      else if (m.type === 'circuit') v.setCircuit(m);
    }

    process(inputs, outputs) {
      const out = outputs[0][0];
      const n = out.length;
      const v = this.voice;
      // pots move smoothly, about 10 ms, so a fast turn does not step the ladder current
      for (let j = 0; j < KNOB_NAMES.length; j++) { const k = KNOB_NAMES[j]; v.knob[k] += (this.target[k] - v.knob[k]) * 0.25; }
      const t0 = clock();
      this.seq.plan(this.frame, n);
      v.process(out, n);
      const dt = clock() - t0;
      this.busyMs += dt;
      if (dt > this.worstMs) this.worstMs = dt;
      if (dt > this.budgetMs) this.late++;
      for (let c = 1; c < outputs[0].length; c++) outputs[0][c].set(out);
      this.frame += n;
      this.peakUa = Math.max(this.peakUa, v.IL * 1e6);
      this.meterCount += n;
      if (this.seq.current !== this.lastStep || this.meterCount >= this.meterEvery) {
        const load = this.busyMs / (1000 * this.meterCount / sampleRate);
        this.port.postMessage({ type: 'state', step: this.seq.current, ua: this.peakUa, load, worstMs: this.worstMs, late: this.late });
        this.worstMs = 0;
        this.late = 0;
        this.lastStep = this.seq.current;
        this.meterCount = 0;
        this.busyMs = 0;
        this.peakUa = 0;
      }
      return true;
    }
  }
  registerProcessor('x0x-vcf', X0xProcessor);
}

if (typeof module !== 'undefined' && module.exports) module.exports = api;
