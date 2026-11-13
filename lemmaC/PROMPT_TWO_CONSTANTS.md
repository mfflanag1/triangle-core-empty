# Launch prompt: two-constants conjecture

Paste the block below into a new session started in
`/Users/maxf/projects/economics-research/game-theory/coopgpt`.

---

Use multi-agent orchestration aggressively and dynamically: fan out subagents, or run a
workflow. You may run up to 8 agents at once and about 40 agent-runs in total, which is
roughly five rounds. Do not use a fixed assignment such as "N agents for strategy X".
Manage the search with the heuristics below.

## Goal

Prove the two-constants conjecture for the random triangle market, or disprove it.

**Setting.** Three sectors of q agents. The pair values A[a,b], B[b,c], C[a,c] are 3q²
independent U[0,1] variables. A team is one agent from each sector and is worth
v(a,b,c) = A[a,b] + B[b,c] + C[a,c]. T is the best total over perfect assignments into q
teams. LP* is the optimum of the fractional relaxation. D = LP* − T. The core of the market
is nonempty exactly when D = 0.

**Conjecture.** There is c > 0 with E D ≥ c·q^(1/3) for all large q. Equivalently, the
per-agent loss constants differ: 3q − LP* ~ a_LP·q^(1/3) and 3q − T ~ a_IP·q^(1/3) with
a_LP < a_IP.

**What is enough.** Anything that gives E D / q^(1/6) → ∞ proves that
P(core nonempty) → 0, through Proposition R: P(core nonempty) ≤ 8·E(3q − T)/(E D)², with
E(3q − T) = Θ(q^(1/3)) (Frieze–Pegden–Tkocz, arXiv:1901.07167, Theorem 2).

## Where things stand

Everything is in `average_case_protection_2026-11-09/`. Read
`lemmaC/SYNTHESIS_2026-11-11.md` first, then the files it names. None of it has been
checked by a different model family, so re-verify any lemma before you build on it.

**Common brief: give this to every agent.**
- The setting, the conjecture and Proposition R, as above.
- Var T and Var LP* are at most 4·E(3q − T) (Poincaré inequality in exponential
  coordinates, X = 1 − exp(−E)).
- With high probability 3q − T ≥ 2.15·q^(1/3). Also E(3q − LP*) ≥ 1.36·q^(1/3).
- D > 0 forces an odd Berge cycle of teams in the support of the LP optimum.
- If only the order of the pair values mattered, the core would never be empty. A proof
  has to use how close the values are to each other, at scale q^(−2/3).
- Conditioning on the optimal assignment lowers LP* (Harris), so treating the unused pair
  values as fresh overestimates LP*.
- Data up to q = 64: E D rises from 0.039 (q = 16) to 0.212 (q = 64) with standard
  deviation about 0.05; the fit is about 0.11·q^(1/3) − 0.25. The LP optimum differs from
  the integer optimum on one giant connected structure that touches about 80% of the
  optimal teams.
- Routes already tried, with the obstacle each one hit:
  - Separating the constants with existing bounds: the best rigorous bounds are 18% apart
    and the true constants differ by about 2%.
  - Improving structures of bounded size: their expected number goes to 0.
  - Summing sub-gaps over blocks of the optimal assignment: blocks of up to half the teams
    capture almost none of D.
  - Pair-price certificates: they exist with probability exactly 1/q!.
  A new idea that revives one of these is welcome if it says why the obstacle no longer
  applies.

**For you, the orchestrator, only.** The currently favoured approach is a local-limit
program:
1. Show that the near-tight teams, at scale q^(−2/3), converge to a Poisson 3-uniform
   hypertree in which each agent's team costs have intensity ξ²/2.
2. Show that on that tree the integral solution is not the cheapest invariant fractional
   perfect matching.
3. Carry the limit constants back to finite q.

A non-rigorous cavity computation gives a_LP ≈ 2.197 and a_IP ≥ 2.2495
(`lemmaC/cavity/CAVITY.md`). A useful first lemma: with high probability the optimal
assignment uses no pair with deficit above q^(−2/3+ε).

## Heuristics

- **Start with a diverse portfolio.** Agents should explore substantially different
  formulations: local weak convergence and the objective method; interpolation and
  comparison arguments; dual certificates and why none can exist; first- and second-moment
  counts over structures whose size grows with q; two-round exposure and sprinkling;
  correlation inequalities; recursive distributional equations; an easier sister model
  (independent team values, or a sparse random 3-partite hypergraph) with a transfer
  argument; exact small cases and the geometry of the assignment polytope; and
  computational sanity checks. Eight agents cannot cover all of these at once, so start
  with the ones that differ most from each other and keep the rest for later rounds.
- **Preserve independence.** Do not tell most agents the favoured approach or the cavity
  numbers. In the first round, at most two of the eight agents work on the local-limit
  program. The rest get the common brief only, so that they do not all converge on the
  same attractive but incomplete reduction.
- **Reallocate after every round.** Read every report. Stop a line that has restated the
  problem without producing a new lemma. Split a line that produced a verified lemma. When
  two independent lines arrive at the same intermediate statement, treat that as a signal
  and give the statement its own agent. Pass verified lemmas between lines; do not pass
  motivations or hunches.
- **Attack every claim.** A lemma claimed as complete goes to a fresh agent that sees only
  the statement and the proof and is told to find the error. Test intermediate claims
  numerically wherever that is possible: `lemmaC/structure/mustruct.py` and
  `lemmaC/scaling/cg.py` compute LP*, T and D; the machine is shared, so use at most 4
  worker processes. A lemma counts as verified only after it survives. Keep about ten of
  the forty runs for these checks. A full proof needs
  two independent adversarial reads, plus a read by Codex through
  `/Users/maxf/.local/bin/cross-model` (do not run its `qa` or `adjudicate` steps).
- **Keep one agent arguing the other side.** The prior estimate that the conjecture is
  true is about 60%. One agent per round tries to show E D = o(q^(1/3)), or finds the
  weakest point of the leading line.
- **Rank partial results.** If the conjecture does not fall, these are worth having, in
  this order:
  1. E D → ∞ at any rate, together with the pair-deficit lemma above. That makes both
     variances o(1), which is enough for the theorem.
  2. limsup P(core nonempty) < 1.
  3. Step 2 of the local-limit program as a standalone statement about the tree.
  4. The conjecture in a sister model.

## Records and rules

- Keep `lemmaC/round3/CLAIMS.md`: every statement with its status (conjectured, claimed,
  verified or refuted), the file that holds it and who checked it. Each agent writes in its
  own subfolder of `lemmaC/round3/`.
- Label a step complete only if you would defend it to a referee. Label everything else as
  a gap. Floating-point results are evidence, never proof. Cite only sources that an agent
  has opened.
- Write only inside the lane folder. Do not commit and do not run `git add`. No paid API
  calls, and no outreach, posting or submission.
- When you finish, append one dated entry to the lane's `LEDGER.md` and record a pending
  receipt with `/Users/maxf/.codex/skills/research-handoff/SKILL.md`.

## When to stop

Stop when a proof survives verification, when two rounds in a row produce no new verified
lemma, or when the agent budget is spent. Then report: what is proved and verified, what
is claimed but unverified, what was refuted, the exact remaining gap, and your confidence
in each.
