#!/usr/bin/env python3
"""Write checks/lean/ErdosSzekeresCertificates.lean from the witness and certificate
files in highdim/.  The Lean file is committed; running this script reproduces it.
Usage: python3 checks/lean/make_lean.py [--stdout]"""
import json
import os
import sys
from itertools import permutations

HERE = os.path.dirname(os.path.abspath(__file__))
HD = os.path.join(HERE, "..", "..", "highdim")


def parity_list(m):
    out = []
    for s in permutations(range(m)):
        inv = sum(1 for i in range(m) for j in range(i + 1, m) if s[i] > s[j])
        out.append(f"([{', '.join(map(str, s))}], {'true' if inv % 2 == 0 else 'false'})")
    return out


def lst(xs):
    return "[" + ", ".join(map(str, xs)) + "]"


def block(name, P, R):
    pts = ",\n  ".join(lst(x) for x in P)
    rels = ",\n  ".join(f"({r['p']}, {lst(r['simplex'])}, {lst(r['weights'])})" for r in R)
    return (f"/-- The points of {name}, indexed from 0 (the paper's $x_{{i+1}}$ is entry $i$). -/\n"
            f"def {name} : List (List Nat) := [\n  {pts}]\n\n"
            f"/-- The relations of {name}: $(p, \\sigma, m)$ with $\\sum_j m_j x_{{\\sigma_j}} = M x_p$. -/\n"
            f"def {name}rels : List (Nat × List Nat × List Nat) := [\n  {rels}]\n")


# (name, witness file, certificate file, d, k, description): the point sets of the paper
SETS = [
    ("X4", "witness_r4_11.json", "cert_r4_11.json", 4, 8, "eleven points in R^4 with no 8 in convex position"),
    ("X3", "witness_r3_15.json", "cert_r3_15.json", 3, 8, "fifteen points in R^3 with no 8 in convex position"),
    ("Y4", "witness_r4_13_k9.json", "cert_r4_13_k9.json", 4, 9, "thirteen points in R^4 with no 9 in convex position"),
    ("Y3", "witness_r3_18_k9.json", "cert_r3_18_k9.json", 3, 9, "eighteen points in R^3 with no 9 in convex position"),
]


def theorems(name, P, R, d, k, desc):
    n = len(P)
    perms = "perms5" if d == 4 else "perms4"
    # the general position check is split into pieces of c simplices each (kernel memory)
    from math import comb
    c = 200 if d == 4 else 400
    m = (comb(n, d + 1) - 1) // c  # number of full pieces before the last one
    L = f"(combs (List.range {name}.length) ({d} + 1))"
    P_ = f"(nondegenerate {perms} {name})"
    pieces = ""
    for i in range(m + 1):
        lst_ = L
        for _ in range(i):
            lst_ = f"({lst_}.drop {c})"
        stmt = f"({lst_}.take {c}).all {P_} = true" if i < m else f"{lst_}.all {P_} = true"
        pieces += f"\ntheorem {name}_gp_{i} : {stmt} := by decide +kernel\n"
    proof = f"{name}_gp_{m}"
    for i in reversed(range(m)):
        proof = f"all_split _ {c} _ {name}_gp_{i}\n    ({proof})"
    return pieces + f"""
theorem {name}_general_position : generalPosition {perms} {name} {d} = true :=
  {proof}

theorem {name}_relations : {name}rels.all (relValid {perms} {name} {d}) = true := by decide +kernel

theorem {name}_size : {name}.length = {n} ∧ {name}rels.length = {len(R)} := by decide +kernel

/-- Every list of at most {n - k} points misses the support of some relation of {name};
equivalently, every {k}-subset of the {desc.split(' points')[0]} points contains such a support. -/
theorem {name}_cover : ∀ V : List Nat, V.length ≤ {n - k} →
    ∃ r, r ∈ {name}rels ∧ ∀ y, y ∈ r.1 :: r.2.1 → y ∉ V := by
  intro V hV
  have ⟨t, ht, htV⟩ := search_spec {n - k} [] ({name}rels.map supp) (by decide +kernel) V hV
    (fun _ _ h => nomatch h)
  have ⟨r, hr, hrt⟩ := mem_map_inv supp t {name}rels ht
  exact ⟨r, hr, fun y hy => htV y (hrt ▸ (show y ∈ supp r from hy))⟩
"""


def main():
    data = [(name, json.load(open(os.path.join(HD, w))), json.load(open(os.path.join(HD, c))), d, k, desc)
            for name, w, c, d, k, desc in SETS]
    listing = "\n".join(f"    {name}: {desc} ({len(R)} relations);"
                        for name, P, R, d, k, desc in data)
    head = r'''/-
  Certificates for "New Lower Bounds for Erdos-Szekeres Numbers in Three and Four
  Dimensions", checked by the Lean 4 kernel.  No Mathlib; core Lean only.

  The point sets:
''' + listing + r'''

  For each set X of n points in R^d with no k points in convex position this file proves:

  * general position: every d+1 points have a nonzero determinant det(1, x);
  * every relation (p, sigma, m) has positive integer weights m_j, satisfies the identity
    sum_j m_j x_{sigma_j} = M x_p with M = sum_j m_j coordinate by coordinate, and has
    det(1, x_{sigma}) /= 0, so x_p lies in the interior of the simplex sigma (Lemma II.3);
  * the covering condition of Proposition II.5: every list V of at most n - k points
    misses the whole support {p} u sigma of some relation.  Taking for V the n - k points
    outside a k-subset T, this says that T contains a support.

  The determinants are computed by the Leibniz formula, split into the sum over even and
  the sum over odd permutations, so that all arithmetic is on natural numbers (all
  coordinates are nonnegative): det /= 0 means that the two sums differ.  The covering
  condition is decided by a branching search, `search`, whose soundness is the theorem
  `search_spec`: the first support either misses V or meets V in a point e, and then V
  without e must miss a support that does not contain e; each tried point e is excluded
  from the later branches.

  Generated by checks/lean/make_lean.py from highdim/*.json.  Check with
    lean checks/lean/ErdosSzekeresCertificates.lean
  The finite checks are proved by `decide +kernel`; `#print axioms` lists what each theorem uses.
-/

/-- Permutations of `0..3` with their parity (`true` for even). -/
def perms4 : List (List Nat × Bool) := [
  ''' + ",\n  ".join(parity_list(4)) + r''']

/-- Permutations of `0..4` with their parity (`true` for even). -/
def perms5 : List (List Nat × Bool) := [
  ''' + ",\n  ".join(parity_list(5)) + r''']

def sumList (l : List Nat) : Nat := l.foldl (· + ·) 0

/-- The even and odd parts of the Leibniz expansion of the determinant of `A`. -/
def detParts (perms : List (List Nat × Bool)) (A : List (List Nat)) : Nat × Nat :=
  perms.foldl (fun acc s =>
    let t := ((List.range s.1.length).zip s.1).foldl
      (fun prod ij => prod * ((A.getD ij.1 []).getD ij.2 0)) 1
    if s.2 then (acc.1 + t, acc.2) else (acc.1, acc.2 + t)) (0, 0)

/-- The matrix with rows `(1, x_s)` for `s` in `S`. -/
def rowsOf (X : List (List Nat)) (S : List Nat) : List (List Nat) :=
  S.map (fun s => 1 :: X.getD s [])

/-- `det(1, x_s : s ∈ S) ≠ 0`. -/
def nondegenerate (perms : List (List Nat × Bool)) (X : List (List Nat)) (S : List Nat) : Bool :=
  let D := detParts perms (rowsOf X S)
  D.1 != D.2

/-- All `k`-element sublists of a list. -/
def combs : List Nat → Nat → List (List Nat)
  | _, 0 => [[]]
  | [], _ + 1 => []
  | x :: xs, k + 1 => (combs xs k).map (x :: ·) ++ combs xs (k + 1)

def generalPosition (perms : List (List Nat × Bool)) (X : List (List Nat)) (d : Nat) : Bool :=
  (combs (List.range X.length) (d + 1)).all (nondegenerate perms X)

/-- A relation `(p, σ, m)`: `d + 1` positive weights, `p ∉ σ`, the identity
`∑ m_j x_{σ_j} = M x_p` in every coordinate, and a nondegenerate simplex `σ`. -/
def relValid (perms : List (List Nat × Bool)) (X : List (List Nat)) (d : Nat)
    (r : Nat × List Nat × List Nat) : Bool :=
  let p := r.1
  let S := r.2.1
  let m := r.2.2
  let M := sumList m
  S.length == d + 1 && m.length == d + 1 && m.all (· > 0) && !S.contains p &&
  (List.range d).all (fun t =>
    sumList ((m.zip S).map (fun ms => ms.1 * (X.getD ms.2 []).getD t 0)) == M * (X.getD p []).getD t 0) &&
  nondegenerate perms X S

/-- `has x l`: does `x` occur in `l`? -/
def has (x : Nat) : List Nat → Bool
  | [] => false
  | y :: ys => Nat.beq y x || has x ys

/-- `rm1 x l`: the list `l` with its first occurrence of `x` removed. -/
def rm1 (x : Nat) : List Nat → List Nat
  | [] => []
  | y :: ys => cond (Nat.beq y x) ys (y :: rm1 x ys)

/-- `minus X s`: the points of `s` that are not in `X`. -/
def minus (X : List Nat) : List Nat → List Nat
  | [] => []
  | y :: ys => cond (has y X) (minus X ys) (y :: minus X ys)

/-- `dropHas e S`: the supports in `S` that do not contain `e`. -/
def dropHas (e : Nat) : List (List Nat) → List (List Nat)
  | [] => []
  | s :: S => cond (has e s) (dropHas e S) (s :: dropHas e S)

/-- Runs `f X (dropHas e S)` for the successive candidates `e`, adding each tried
candidate to the excluded points `X`. -/
def branchWith (f : List Nat → List (List Nat) → Bool) (S : List (List Nat)) :
    List Nat → List Nat → Bool
  | _, [] => true
  | X, e :: es => f X (dropHas e S) && branchWith f S (e :: X) es

/-- `search b X S`: every list `V` of at most `b` points, none of them in `X`, misses
some support in `S` entirely.  The first support `s` of `S` either misses `V`, or `V`
contains a point `e` of `s` outside `X`, and removing `e` from `V` leaves at most `b - 1`
points that must miss a support of `S` not containing `e`. -/
def search : Nat → List Nat → List (List Nat) → Bool
  | 0, _, [] => false
  | 0, _, _ :: _ => true
  | _ + 1, _, [] => false
  | b + 1, X, s :: S => branchWith (search b) (s :: S) X (minus X s)

theorem and_split : ∀ a b : Bool, (a && b) = true → a = true ∧ b = true
  | true, true, _ => ⟨rfl, rfl⟩
  | false, _, h => Bool.noConfusion h
  | true, false, h => Bool.noConfusion h

theorem mem_of_has (x : Nat) : ∀ l : List Nat, has x l = true → x ∈ l
  | [], h => Bool.noConfusion h
  | y :: ys, h => by
    have h1 : (Nat.beq y x || has x ys) = true := h
    revert h1
    cases hy : Nat.beq y x with
    | true => intro _; rw [← Nat.eq_of_beq_eq_true hy]; exact List.Mem.head ys
    | false => intro h1; exact List.Mem.tail y (mem_of_has x ys h1)

theorem not_mem_of_has (x : Nat) : ∀ l : List Nat, has x l = false → x ∉ l
  | [], _, hm => nomatch hm
  | y :: ys, h, hm => by
    have h1 : (Nat.beq y x || has x ys) = false := h
    revert h1
    cases hy : Nat.beq y x with
    | true => intro h1; exact Bool.noConfusion h1
    | false =>
      intro h1
      cases hm with
      | head => rw [Nat.beq_refl] at hy; exact Bool.noConfusion hy
      | tail _ hm' => exact not_mem_of_has x ys h1 hm'

theorem length_rm1 (x : Nat) : ∀ l : List Nat, x ∈ l → (rm1 x l).length + 1 = l.length
  | [], hm => nomatch hm
  | y :: ys, hm => by
    show (cond (Nat.beq y x) ys (y :: rm1 x ys)).length + 1 = ys.length + 1
    cases hy : Nat.beq y x with
    | true => rfl
    | false =>
      cases hm with
      | head => rw [Nat.beq_refl] at hy; exact Bool.noConfusion hy
      | tail _ hm' => exact congrArg (· + 1) (length_rm1 x ys hm')

theorem mem_rm1 (x y : Nat) : ∀ l : List Nat, y ∈ l → y ≠ x → y ∈ rm1 x l
  | [], hm, _ => nomatch hm
  | z :: zs, hm, hne => by
    show y ∈ cond (Nat.beq z x) zs (z :: rm1 x zs)
    cases hz : Nat.beq z x with
    | true =>
      cases hm with
      | head => exact absurd (Nat.eq_of_beq_eq_true hz) hne
      | tail _ hm' => exact hm'
    | false =>
      cases hm with
      | head => exact List.Mem.head _
      | tail _ hm' => exact List.Mem.tail z (mem_rm1 x y zs hm' hne)

theorem mem_of_mem_rm1 (x y : Nat) : ∀ l : List Nat, y ∈ rm1 x l → y ∈ l
  | [], hm => nomatch hm
  | z :: zs, hm => by
    have h1 : y ∈ cond (Nat.beq z x) zs (z :: rm1 x zs) := hm
    revert h1
    cases Nat.beq z x with
    | true => intro h1; exact List.Mem.tail z h1
    | false =>
      intro h1
      cases h1 with
      | head => exact List.Mem.head _
      | tail _ h' => exact List.Mem.tail z (mem_of_mem_rm1 x y zs h')

theorem mem_minus (X : List Nat) (y : Nat) : ∀ s : List Nat, y ∈ s → y ∉ X → y ∈ minus X s
  | [], hm, _ => nomatch hm
  | z :: zs, hm, hX => by
    show y ∈ cond (has z X) (minus X zs) (z :: minus X zs)
    cases hz : has z X with
    | true =>
      cases hm with
      | head => exact absurd (mem_of_has y X hz) hX
      | tail _ hm' => exact mem_minus X y zs hm' hX
    | false =>
      cases hm with
      | head => exact List.Mem.head _
      | tail _ hm' => exact List.Mem.tail z (mem_minus X y zs hm' hX)

theorem mem_dropHas (e : Nat) (t : List Nat) :
    ∀ S : List (List Nat), t ∈ dropHas e S → t ∈ S ∧ has e t = false
  | [], hm => nomatch hm
  | s :: S, hm => by
    have h1 : t ∈ cond (has e s) (dropHas e S) (s :: dropHas e S) := hm
    revert h1
    cases hs : has e s with
    | true =>
      intro h1
      have ⟨h2, h3⟩ := mem_dropHas e t S h1
      exact ⟨List.Mem.tail s h2, h3⟩
    | false =>
      intro h1
      cases h1 with
      | head => exact ⟨List.Mem.head _, hs⟩
      | tail _ h' =>
        have ⟨h2, h3⟩ := mem_dropHas e t S h'
        exact ⟨List.Mem.tail s h2, h3⟩

theorem branch_spec (f : List Nat → List (List Nat) → Bool) (S : List (List Nat)) (b : Nat)
    (hf : ∀ X S', f X S' = true → ∀ V : List Nat, V.length ≤ b → (∀ y, y ∈ V → y ∉ X) →
      ∃ t, t ∈ S' ∧ ∀ y, y ∈ t → y ∉ V) :
    ∀ es X : List Nat, branchWith f S X es = true → ∀ V : List Nat, V.length ≤ b + 1 →
      (∀ y, y ∈ V → y ∉ X) → (∀ y, y ∈ es → y ∉ V) ∨ ∃ t, t ∈ S ∧ ∀ y, y ∈ t → y ∉ V
  | [], _, _, _, _, _ => Or.inl (fun _ hy => nomatch hy)
  | e :: es, X, h, V, hV, hX => by
    have h1 : (f X (dropHas e S) && branchWith f S (e :: X) es) = true := h
    have ⟨ha, hb⟩ := and_split _ _ h1
    cases he : has e V with
    | true =>
      have heV := mem_of_has e V he
      have hlen : (rm1 e V).length ≤ b := by
        rw [← length_rm1 e V heV] at hV
        exact Nat.le_of_succ_le_succ hV
      have hX' : ∀ y, y ∈ rm1 e V → y ∉ X := fun y hy => hX y (mem_of_mem_rm1 e y V hy)
      have ⟨t, htS, htV⟩ := hf X (dropHas e S) ha (rm1 e V) hlen hX'
      have ⟨htS', hte⟩ := mem_dropHas e t S htS
      refine Or.inr ⟨t, htS', fun y hy hyV => ?_⟩
      cases Nat.decEq y e with
      | isTrue hye => rw [hye] at hy; exact not_mem_of_has e t hte hy
      | isFalse hye => exact htV y hy (mem_rm1 e y V hyV hye)
    | false =>
      have heV := not_mem_of_has e V he
      have hX' : ∀ y, y ∈ V → y ∉ e :: X := fun y hy hyeX => by
        cases hyeX with
        | head => exact heV hy
        | tail _ h' => exact hX y hy h'
      cases branch_spec f S b hf es (e :: X) hb V hV hX' with
      | inl hl =>
        refine Or.inl (fun y hy => ?_)
        cases hy with
        | head => exact heV
        | tail _ h' => exact hl y h'
      | inr hr => exact Or.inr hr

theorem search_spec : ∀ (b : Nat) (X : List Nat) (S : List (List Nat)), search b X S = true →
    ∀ V : List Nat, V.length ≤ b → (∀ y, y ∈ V → y ∉ X) → ∃ t, t ∈ S ∧ ∀ y, y ∈ t → y ∉ V
  | 0, _, [], h, _, _, _ => Bool.noConfusion h
  | 0, _, s :: _, _, V, hV, _ => ⟨s, List.Mem.head _, fun _ _ hyV => by
      cases V with
      | nil => exact nomatch hyV
      | cons _ _ => exact absurd hV (Nat.not_succ_le_zero _)⟩
  | _ + 1, _, [], h, _, _, _ => Bool.noConfusion h
  | b + 1, X, s :: S, h, V, hV, hX => by
    cases branch_spec (search b) (s :: S) b (search_spec b) (minus X s) X h V hV hX with
    | inl hl => exact ⟨s, List.Mem.head _, fun y hy hyV => hl y (mem_minus X y s hy (hX y hyV)) hyV⟩
    | inr hr => exact hr

theorem mem_map_inv {α β : Type} (f : α → β) (b : β) : ∀ l : List α, b ∈ l.map f → ∃ a, a ∈ l ∧ f a = b
  | [], h => nomatch h
  | a :: l, h => by
    cases h with
    | head => exact ⟨a, List.Mem.head _, rfl⟩
    | tail _ h' =>
      have ⟨a', ha', hfa⟩ := mem_map_inv f b l h'
      exact ⟨a', List.Mem.tail a ha', hfa⟩

/-- The support `{p} ∪ σ` of a relation `(p, σ, m)`. -/
def supp (r : Nat × List Nat × List Nat) : List Nat := r.1 :: r.2.1

/-- A list passes `all f` if its first `n` entries and the rest do; the general position
checks are split in this way into pieces that the kernel checks one at a time. -/
theorem all_split {α : Type} (f : α → Bool) : ∀ (n : Nat) (l : List α),
    (l.take n).all f = true → (l.drop n).all f = true → l.all f = true
  | 0, _, _, h => h
  | _ + 1, [], _, _ => rfl
  | n + 1, a :: l, h1, h2 => by
    have ⟨ha, ht⟩ := and_split (f a) ((l.take n).all f) h1
    have hl := all_split f n l ht h2
    show (f a && l.all f) = true
    rw [ha, hl]
    rfl

'''
    body = "\n".join(block(name, P, R) for name, P, R, d, k, desc in data)
    thms = "".join(theorems(*x) for x in data)
    axioms = "\n".join(f"#print axioms {name}_{t}" for name, *_ in data
                       for t in ("general_position", "relations", "size", "cover"))
    text = head + body + thms + "\n" + axioms + "\n"
    if "--stdout" in sys.argv:
        sys.stdout.write(text)
    else:
        open(os.path.join(HERE, "ErdosSzekeresCertificates.lean"), "w").write(text)


if __name__ == "__main__":
    main()
