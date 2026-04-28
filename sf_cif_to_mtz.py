#!/usr/bin/env python3

import argparse
from pathlib import Path

import gemmi


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="Convert a reflection block from an SF CIF file to MTZ."
	)
	parser.add_argument("input_cif", help="Input structure-factor CIF file path")
	parser.add_argument(
		"-o",
		"--output",
		help="Output MTZ file path (default: same name as input with .mtz extension)",
	)
	parser.add_argument(
		"-b",
		"--block-index",
		type=int,
		default=0,
		help="Reflection block index to extract (default: 0)",
	)
	return parser.parse_args()


def main() -> None:
	args = parse_args()

	input_path = Path(args.input_cif)
	output_path = Path(args.output) if args.output else input_path.with_suffix(".mtz")

	doc = gemmi.cif.read_file(str(input_path))
	refl_blocks = gemmi.as_refln_blocks(doc)

	if not refl_blocks:
		raise SystemExit(f"No reflection blocks found in {input_path}")
	if args.block_index < 0 or args.block_index >= len(refl_blocks):
		raise SystemExit(
			f"Block index {args.block_index} out of range (available: 0..{len(refl_blocks)-1})"
		)

	selected_block = refl_blocks[args.block_index]

	# Gemmi API compatibility: newer versions expose ReflnBlock.make_mtz(),
	# while older versions use CifToMtz.convert_block_to_mtz().
	mtz = gemmi.CifToMtz().convert_block_to_mtz(selected_block)

	mtz.write_to_file(str(output_path))
	print(f"Wrote {output_path}")


if __name__ == "__main__":
	main()