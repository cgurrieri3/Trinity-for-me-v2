# M3Long — handoff

**Status: FIXED AND APPLIED 2026-09-01.** Everything below described the broken state and the
proposed fix; it is kept as the rationale. What actually shipped, and the one place this document
was wrong, is summarised here first.

## The source this analysis follows

The reference is a MAGIC thesis chapter (6.3, Head-tail Discrimination), which defines **two**
variables. Reading the old code against it changes the diagnosis below in two places:

    M3Long     = ( (1/N) * sum_i q_i   (x_i - x_cog)^3 )^(1/3)      Eq 6.2
    M3LongPow3 = ( (1/N) * sum_i q_i^3 (x_i - x_cog)^3 )^(1/3)
    "The x axis is thereby the major axis of the Hillas ellipse."

The old code computed `(1/N) sum q^3 dx^3` along **camera x**. That is `M3LongPow3`, not a
corrupted `M3Long`. So of complaint #3's three sub-points below:

- *"Charge is cubed"* -- **not a bug.** That is `M3LongPow3`, chosen deliberately for low-energy
  roundish images where the faint pixels carrying the tail are noise-dominated.
- *"Normalised by pixel count, not charge"* -- **not a bug.** Eq 6.2 normalises by N. And it
  cannot matter for head-tail anyway: both denominators are positive and cbrt preserves sign, so
  1/N vs 1/sum(q) rescales the value and can never change its sign.
- *"No projection onto the major axis"* -- **genuinely wrong.** The thesis is explicit that x is
  the major axis; the old code took separate moments along camera x and camera y. Together with
  the copy-paste and NaN bugs, it implemented neither thesis variable.

## Applied 2026-09-01

- `getM3Long` rewritten as a proper longitudinal third moment: takes `anglerad`, projects onto the
  major axis, returns a `double` via `std::cbrt`. Debug `cout` calls removed.
- **Both thesis variables are computed and stored**, sharing one `projectOnMajorAxis` helper:
  `getM3Long` (q^1, normalised by sum q) and `getM3LongPow3` (q^3, normalised by N). Pixels with
  q <= 0 are dropped in both -- in the q^1 sum a negative charge flips the weighting, in the q^3
  sum it survives the cube with the wrong sign. A perfectly vertical major axis returns 0 rather
  than a coin flip.
- `EventInfo::M3Longx`/`M3Longy` (two `std::string`s) replaced by `Float_t M3Long` and
  `Float_t M3LongPow3`. `ClassDef` bumped 2 -> 4. Both copies of `EventInfo.{h,cpp}`
  (ClusterCleaning and EventSelection) updated together, dictionaries regenerated, both programs
  rebuilt clean from scratch.
- Canvas labels in `EventCleaning.h` and `eventSelection.cpp` print `M3Long:%.2f M3Pow3:%.3g`.
- Verified: `Claude/test_m3long.cpp` (synthetic blobs, exact mirror antisymmetry, vanishing moment
  for a symmetric blob, no NaN on degenerate input), and a real night (20250517): 0 NaN anywhere in
  the 25867-entry tree, and 6 positive / 4 negative `M3Long` over the 10 surviving events rather
  than 100% "+".

### The two variables carry OPPOSITE signs -- do not mix them

This is the practical trap. Cubing the charge lets the bright compact head dominate the sum
instead of the faint extended tail, so `M3LongPow3` generally reports the opposite direction to
`M3Long` on the same image. Measured, not assumed:

- Synthetic blobs: opposed at every angle tested (0, +/-30, -45, 60, 89 deg).
- Real night 20250517: **7 of the 10 surviving events opposed**, 3 agreed.

So a head-tail cut must be written against one variable and stay with it. The sign discussion in
the next section applies to `M3Long` only.

### Correction: the sign does NOT point to the brighter end

The sections below claim `M3Long > 0` means the brighter end is toward +x. That is wrong, and the
formula this document prescribes does not do it. Because `(xcog,ycog)` is the *charge-weighted*
centroid, piling charge at one end drags the centroid to that end and leaves a long faint tail at
the other, and the cubed distance makes that tail dominate the sum. The sign is therefore the
direction of the **skew**, pointing away from the bright end.

Hand-checked for a blob at 0 deg with `d = -4..4` and `q = 90 + 20d` (charge piled toward +d):
centroid at `d = 1.4815`, `sum q (d-dbar)^3 / sum q = -5.6451`, `cbrt = -1.7806`, which is exactly
what the shipped code returns. So `M3Long > 0` means the image tails off toward +x and the denser,
brighter end sits toward **-x**.

This is the standard Hillas convention, and matches this document's own complaint #3 that the old
code "is not the standard quantity" — so the formula was kept and the interpretation corrected. If
what you actually want is a quantity whose sign points at the bright end, negate it at the point of
use rather than inside `getM3Long`.

---

*Original document follows, written 2026-08-31.*

## Why this came up

We wanted to know which end of the shower ellipse holds more signal. The major-axis angle cannot
answer that: PCA is built from *second* moments, which are symmetric about the COG, and
`atan(eigenVecs[1][0]/eigenVecs[0][0])` returns (-90°, +90°), so it describes a headless line. The
sign of the angle only means "tilts up to the right" vs "down to the right".

That asymmetry is exactly what a *third* moment is for, which is what `M3Long` is supposed to
provide. It currently does not.

(The angle itself is fine. Verified: a synthetic blob built at +30° returns `30.00°`, and ROOT's
`TMatrixDSymEigen` sorts eigenvalues descending, so column 0 really is the major axis and `L`/`W`
are assigned correctly.)

## Where the code is

| What | Location |
|---|---|
| Declaration | `EventCleaning.h:162` |
| Definition | `EventCleaning.h:582-622` |
| Call site | `EventCleaning.h:501` |
| Sign classification | `EventCleaning.h:502-513` |
| Stored on EventInfo | `SetM3Longx` / `SetM3Longy`, as `std::string` (`EventInfo.h:130-131`) |
| Consumed | PDF label only — `EventCleaning.h:515-522`, `eventSelection.cpp:357-364` |

It is never plotted and never cut on. Nothing depends on the numeric value, only on the `"+"`/`"-"`/`"0"`
string in a canvas title. **That makes this cheap to change.**

## What is wrong

### 1. Copy-paste error — the y component is a function of x

`EventCleaning.h:616-617`:

```cpp
M3Long[0] = std::pow((M3Long[0]/psize),0.333333333);
M3Long[1] = std::pow((M3Long[0]/psize),0.333333333);  // M3Long[0], should be [1]
```

`M3Long[0]` is overwritten on the first line, so the second line consumes the already-reduced x
value. `M3Longy` is meaningless.

### 2. Every negative skew is reported as `"+"`

`std::pow` with a negative base and a fractional exponent returns NaN. Confirmed in the container:

```
pow(-8, 0.333333333) = -nan
```

Skewness is signed — negative values are the entire point of the quantity. And the classifier at
`EventCleaning.h:502-513` tests `< 0` then `== 0`, both of which are **false** for NaN, so it falls
through to the `"+"` initialiser. Every negative-skew event is labelled `"+"`.

So the one field meant to encode "which way" reports `"+"` for roughly half the events by
construction, silently.

### 3. It is not the standard quantity

Three deviations from the Hillas third moment, all in `EventCleaning.h:598-616`:

- **Charge is cubed.** `std::pow(amps[...],3.0) * std::pow(dx,3.0)` — the standard weights by `q¹`,
  giving `Σ q·d³`. Cubing the charge makes the brightest pixel dominate to the ninth power overall.
- **No projection onto the major axis.** It takes separate third moments along camera x and y. The
  "Long" in M3Long means *longitudinal* — along the major axis. As written it is not that.
- **Normalised by pixel count, not charge.** `M3Long[0]/psize` should be `/ Σ q_i`.

## What a fix looks like

Project each surviving pixel onto the major-axis unit vector, weight by charge, normalise by total
charge, and take a sign-preserving cube root:

```
d_i     = (x_i - xcog)·cosθ + (y_i - ycog)·sinθ
M3Long  = cbrt( Σ q_i · d_i³ / Σ q_i )
```

```cpp
// θ is anglerad, already computed in CompletePanel4 — getM3Long needs it passed in.
double getM3Long(double xcog, double ycog, double anglerad,
                 std::vector<int> sur_pix, std::vector<float> amps) {
    double num = 0.0, qsum = 0.0;
    double ct = cos(anglerad), st = sin(anglerad);
    for (std::vector<int>::size_type p = 0; p < sur_pix.size(); p++) {
        int nx, ny;
        plottools->FindBin(sur_pix[p], &nx, &ny);
        double q = amps[sur_pix[p]];
        if (q <= 0) continue;                       // negative charge would flip the weighting
        double d = (nx - xcog)*ct + (ny - ycog)*st; // signed distance along the major axis
        num  += q * d*d*d;
        qsum += q;
    }
    if (qsum <= 0) return 0.0;                      // no charge, no moment
    return std::cbrt(num / qsum);                   // cbrt, NOT pow — keeps the sign
}
```

`std::cbrt(-8) = -2`, which is the whole reason to use it over `std::pow`.

### The sign convention is well defined

Worth stating explicitly, because it is the part most likely to cause confusion: `atan` returns
(-90°, +90°), so `cos θ > 0` always, and the major-axis unit vector always points toward **+x**.
So `M3Long > 0` means the brighter end is toward +x (right side of the camera), `M3Long < 0` means
toward -x. That is consistent across every event.

Only exception: a perfectly vertical axis (`eigenVecs[0][0] == 0`) gives `cos θ = 0`, and the
direction then rides on the arbitrary sign of the eigenvector. Rare, but it should probably be
detected and reported as 0 rather than a coin flip.

## Decisions for whoever picks this up

1. **Replace `M3Longx`/`M3Longy`, or add a new field?** They are `std::string` on `EventInfo` and
   only ever feed a canvas label. Replacing them with a single `Float_t M3Long` is cleaner and
   makes the value plottable, but it is an `EventInfo` schema change — `ClassDef` was already
   bumped to 2 on 2026-08-31 for `SaturatedPixelIDs`, so this would be 3.
2. **Keep an x/y decomposition at all?** Once it is projected onto the major axis, "longitudinal"
   is a single number. A transverse third moment (`M3Trans`, along the minor axis) is a separate,
   also-meaningful quantity if wanted — but it is not what `M3Longy` currently is.
3. **Existing files.** No stored `M3Longx`/`M3Longy` value is trustworthy, so nothing is lost by
   changing the meaning. Anything re-derived from them needs recleaning regardless.

## How to verify a fix

- **Synthetic test.** Build a blob with charge deliberately piled at one end of a known axis;
  confirm the sign points to the heavy end, and that mirroring the blob flips the sign.
- **NaN count must be zero.** The current implementation's failure mode is silent, so assert on it.
- **Sign balance.** Over a real night, a working M3Long should not be ~100% `"+"`. That ratio is the
  fastest smoke test that the NaN bug is gone.
- The debug `cout` calls at `EventCleaning.h:600, 605, 608-609` fire once per pixel per event (plus
  `591` once per event). They should come out while you are in there.

## Related

`Distance2MajorAxis` on `EventInfo` is a single value overwritten per pixel by the cleaning, so it
holds only the last pixel's distance. `eventSelection.cpp` already works around this by recomputing
per-pixel distances itself. Same neighbourhood of code, same kind of problem — worth looking at in
the same pass.
