# Project proposal

## 1. The problem and goal

**The problem and goal in my project:**

A particle in a one-dimensional potential V(x) has a heat capacity Cv. It can be calculated quantum mechanically (Cv_q) or classically (Cv_cl). The two are usually expected to satisfy Cv_q < Cv_cl, but some potentials show Cv_q > Cv_cl in a window of temperatures. This matters for quantum heat engines, because the work an engine produces depends on Cv.

The excess I know of comes from one structure: a narrow trap with few, widely spaced levels inside a wide reservoir with many closely spaced levels. An example is the well −V0/cosh²(αx) inside a box. My working claim, currently a proof sketch and not a finished proof, makes two statements:

- **Sufficiency:** a trap plus a large enough reservoir gives an excess, in a temperature window that can be predicted.
- **Necessity (conjecture):** a smooth potential without this structure gives no excess.

**The goal** is to search for counterexamples to this claim:

| Counterexample | Meaning |
|---|---|
| A smooth potential classified as having no trap plus reservoir, with a confirmed excess above the steep-wall effect | A new family of potentials, which breaks necessity |
| A trap plus reservoir where the predicted window is missing or wrong | The sufficiency criterion must be corrected |

The steep-wall effect is a separate, known excess that needs no trap. It reaches k_B/16 for an empty box, and my tests with steep walls approach that value from below. For now, a result only counts as a new family if it is larger than k_B/16. This is a practical threshold, not a proven bound for all potentials, and it would miss a real effect smaller than k_B/16.

**How I will measure improvement:**

- the number of confirmed counterexamples found, including zero;
- on a task with a known answer, how often the search finds the excess, and how many simulator calls it needs;
- how those numbers compare with random search and a numerical optimizer given the same number of simulator calls;
- on a task with no expected excess, how often the search correctly reports none.

**Limits and requirements:**

- The potentials must be one-dimensional and bounded from below, so V(x) has a finite minimum and no singular wells such as delta functions or −1/x².
- The potentials must also be smooth and confining, within a fixed energy scale and a fixed ratio between the widest and narrowest features. Without these limits the problem is too easy, because the excess keeps growing as the box gets larger. The values are still to be chosen.
- Only numbers calculated by the simulator count as results. The language model never calculates Cv itself.
- A numerical error must not count as a discovery, so every candidate must pass the checks in section 2.
- Finding no counterexample supports the claim but does not prove it.

## 2. The simulator

**The simulator in my project:**

I will write a small Python simulator based on my existing research code. It runs on my laptop with no GUI, is deterministic, and accepts a potential and its settings as a JSON file.

**Inputs:** the potential, built from a list of allowed terms, plus the temperature range and numerical settings.

**The simulator has two tools:**

1. **Screen (called often):**
   - Cv_q comes from a sine-DVR calculation (Colbert and Miller) of the energy levels.
   - Cv_cl comes from the classical integral, where the momentum integral is done by hand: Z_cl = (1/h) ∫∫ exp(−βH) dx dp = (1/λ_th) ∫ exp(−βV(x)) dx, with λ_th = h/√(2πm k_B T). This gives Cv_cl/k_B = 1/2 + β²·Var(V).
   - It also reports the same Cv_q result on a finer grid.
2. **Confirm (finalists only):** the ξ ladder, the same potential solved with ħ/ξ for growing ξ. It reports the value at each ξ and the gap to the integral. It only reports and never stops a run. The ξ values and the agreement tolerance are not set yet, because they depend on the open question in section 3.

**Outputs:** Cv_q, Cv_cl, and ΔCv = Cv_q − Cv_cl at each temperature, the size of the numerical checks, and the quantities used to classify the potential, such as the trap and reservoir weights.

**Speed:** in one test, on a well in a box, the screen step took about 1.3 s for 251 temperatures. Sharp traps need finer grids and cost more, so I will need a limit on grid size. I have not measured that case yet.

**Tests with known answers:** the harmonic oscillator, and −V0/cosh²(αx) in a box, which has closed forms for both Z_q and Z_cl.

## 3. The agent's role

**The agent in my project:**

**What the agent decides:**

- which combination of terms to try next, for example power-law walls, 1/cosh² wells, Gaussians, Woods-Saxon wells, Morse wells and periodic terms;
- the parameter values, within the allowed ranges;
- which claim a candidate is testing, based on what earlier candidates showed;
- whether a surprising result looks real or numerical, and which check to run next;
- which finalists to propose for the expensive confirm step. I approve these at a checkpoint.

The allowed ranges are not chosen yet. They will come from the limits in section 1.

**How it uses the results:** it proposes a potential, the simulator runs it, the agent reads ΔCv and the classification, and it uses them to choose the next potential. For example, if an excess disappears when a term is removed, it can test whether that term was acting as a hidden trap or reservoir.

**What I expect the language model to contribute:** changing the *form* of the potential and working out why a candidate failed, based on the physics in the claim.

**Could a simple program do the same?** Possibly, and that is the main risk. Inside one fixed formula, a numerical optimizer will probably match or beat the agent. So the comparison with random search and an optimizer is part of the evaluation, and I will run it early. If the agent does not beat them, then by the course's definition the project is a script, and I must change the agent's role or the project before the Week 6 revision deadline.

**The one question I still need to answer:** is the classical integral exact for the potentials I will search, so that it can serve as Cv_cl? The momentum integral is exact by hand, but three things still need checking:

- that the integral equals the ħ → 0 limit of the quantum calculation for these potentials, including whether the ξ ladder really approaches it;
- that the remaining x integral is accurate numerically, by comparing integration methods that work differently (adaptive quadrature, fine fixed-grid Simpson, and Gauss-Legendre panels placed around each narrow feature), because a narrow trap can be missed by a routine that never samples it;
- that it matches the known closed forms.

If it is exact, the ξ ladder is only an independent sanity check, and its tolerance does not need to be small.

**Other open points:** a fixed rule for "has no trap plus reservoir" on a general potential, and how many simulator calls and how much API cost each run may use.

## Explain your project to someone outside your field

**My project in simple words:**

Heat capacity tells you how much energy something absorbs when it warms up by one degree. For very small particles there are two ways to calculate it: the "everyday" physics way and the quantum way, where the particle can only hold certain fixed amounts of energy, like standing on the rungs of a ladder. Normally the quantum answer is the smaller one, but for a few special shapes of energy landscape it is larger, which could help tiny engines produce more work. The only shape I know that does this is a small, narrow pit inside a very wide, flat basin. I have a partly proven idea that this is the only kind of shape that works, and I want to try to prove it wrong.

To do that, I have a calculator program that takes a landscape shape and reliably computes both answers, giving the same result every time. The agent is an AI model that acts like a researcher running experiments. It suggests a shape, runs the calculator, reads the result, and decides what to try next. Before I trust any surprising result, a slower and separate calculation double-checks it. I will also test whether the AI does better than a simple program that tries shapes at random. If it does not, the AI is not adding anything.
