from pathlib import Path

from PIL import Image


base_dir = Path("fuzzyevents_1")
source_dir = base_dir / "fuzzy"
original_dir = base_dir / "original_pngs"
flipped_dir = base_dir / "flipped_pngs"

original_dir.mkdir(parents=True, exist_ok=True)
flipped_dir.mkdir(parents=True, exist_ok=True)

for source_path in sorted(source_dir.glob("*.png")):
    original_path = original_dir / source_path.name
    flipped_path = flipped_dir / source_path.name

    with Image.open(source_path) as image:
        image.save(original_path)

        # Flip the already-rendered fuzzy-event PNG in both image axes.
        # This mirrors the same "left <-> right" and "top <-> bottom" change
        # used for the ROOT camera histograms, while preserving Sofia's
        # original PNG in original_pngs.
        image.transpose(Image.ROTATE_180).save(flipped_path)

print(f"Saved originals to {original_dir}")
print(f"Saved flipped PNGs to {flipped_dir}")
